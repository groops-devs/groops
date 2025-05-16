/***********************************************/
/**
* @file slrParametrizationStaticPositions.cpp
*
* @brief Position estimation with no-net constraints.
* @see SlrParametrization
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "inputOutput/system.h"
#include "files/fileInstrument.h"
#include "files/fileGriddedData.h"
#include "misc/varianceComponentEstimation.h"
#include "classes/platformSelector/platformSelector.h"
#include "slr/slrParametrization/slrParametrization.h"
#include "slr/slrParametrization/slrParametrizationStaticPositions.h"

/***********************************************/

SlrParametrizationStaticPositions::SlrParametrizationStaticPositions(Config &config)
{
  try
  {
    sigmaNoNetTranslation = sigmaNoNetRotation = sigmaNoNetScale = 0;
    huber = 1e99; huberPower = 0;

    readConfig(config, "name",                      name,                   Config::OPTIONAL, "parameter.staticPositions", "used for parameter selection");
    readConfig(config, "selectStations",            selectorStations,       Config::MUSTSET,  "", "");
    readConfig(config, "outputfileGriddedPosition", fileNameGrid,           Config::OPTIONAL, "output/gridPosition_{loopTime:%D}.dat",              "delta north east up for all stations");
    readConfig(config, "outputfilePosition",        fileNamePosition,       Config::OPTIONAL, "output/stationPosition_{loopTime:%D}.{station}.dat", "variable {station} available, full estimated coordinates (in TRF)");
    readConfig(config, "nameConstraint",            nameConstraint,         Config::OPTIONAL, "constraint.staticPositions", "used for parameter selection");
    readConfig(config, "selectNoNetStations",       selectorNoNetStations,  Config::DEFAULT,  R"(["all"])", "");
    readConfig(config, "inputfileNoNetPositions",   fileNameNoNetPositions, Config::OPTIONAL, "{groopsDataDir}/slr/stations/position/slrf2020_CM/stationPosition.{station}.dat", "variable {station} available, precise coordinates used for no-net constraints (in TRF)");
    readConfig(config, "noNetTranslationSigma",     sigmaNoNetTranslation,  Config::OPTIONAL, "0",     "(0 = unconstrained) sigma [m] for no-net translation constraint on station coordinates");
    readConfig(config, "noNetRotationSigma",        sigmaNoNetRotation,     Config::OPTIONAL, "0.001", "(0 = unconstrained) sigma [m] at Earth's surface for no-net rotation constraint on station coordinates");
    readConfig(config, "noNetScaleSigma",           sigmaNoNetScale,        Config::OPTIONAL, "0",     "(0 = unconstrained) sigma [m] for no-net scale constraint on station coordinates");
    readConfig(config, "huber",                     huber,                  Config::OPTIONAL, "2.5",   "stations > huber*sigma0 are downweighted in no-net constraint");
    readConfig(config, "huberPower",                huberPower,             Config::OPTIONAL, "1.5",   "stations > huber: sigma=(e/huber)^huberPower*sigma0");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SlrParametrizationStaticPositions::init(Slr *slr, const std::vector<const SlrParametrizationGravityField*> &/*paramGravityField*/)
{
  try
  {
    this->slr  = slr;
    index.resize(slr->stations.size());
    selectedStations = slr->selectStations(selectorStations);

    // apriori positions
    pos.resize(slr->stations.size());
    pos0.resize(slr->stations.size());
    for(auto stat : slr->stations)
      if(selectedStations.at(stat->idStat()) && stat->useable())
        pos.at(stat->idStat()) = pos0.at(stat->idStat()) = stat->pos;

    // no net-positions
    // ----------------
    if(sigmaNoNetRotation || sigmaNoNetTranslation || sigmaNoNetScale)
    {
      std::vector<const Platform*> platforms(slr->stations.size(), nullptr);
      for(auto stat : slr->stations)
        if(selectedStations.at(stat->idStat()) && stat->useable())
          platforms.at(stat->idStat()) = &stat->platform;

      // no net positions
      noNetPos = pos;
      if(!fileNameNoNetPositions.empty())
      {
        noNetPos = std::vector<Vector3d>(pos0.size(), Vector3d(NAN_EXPR, NAN_EXPR, NAN_EXPR));
        Bool someDisabled = FALSE;
        do
        {
          someDisabled = FALSE;
          selectedNoNetStations = selectorNoNetStations->select(slr->times.front(), slr->times.back(), platforms);
          for(UInt idStat=0; idStat<slr->stations.size(); idStat++)
            if(selectedNoNetStations.at(idStat) && std::isnan(noNetPos.at(idStat).r()))
            {
              try
              {
                VariableList fileNameVariableList;
                fileNameVariableList.setVariable("station", slr->stations.at(idStat)->name());
                if(!System::exists(fileNameNoNetPositions(fileNameVariableList)))
                  throw(Exception("file <"+fileNameNoNetPositions(fileNameVariableList).str()+"> not exist"));
                const Time timesMid = 0.5*(slr->times.front()+slr->times.back());
                Vector3dArc arc = InstrumentFile::read(fileNameNoNetPositions(fileNameVariableList));
                auto iter = std::min_element(arc.begin(), arc.end(), [&](const Epoch &e1, const Epoch &e2)
                                            {return std::fabs((e1.time-timesMid).mjd()) < std::fabs((e2.time-timesMid).mjd());});
                if(!arc.size() || ((arc.size() > 1) && (std::fabs((iter->time-timesMid).mjd()) > 0.5*medianSampling(arc.times()).mjd())))
                  throw(Exception("No a-priori position found"));
                noNetPos.at(idStat) = iter->vector3d;
              }
              catch(std::exception &e)
              {
                logWarningOnce<<slr->stations.at(idStat)->name()<<": "<<e.what()<<" -> not used for computation of net translation/rotation/scale"<<Log::endl;
                platforms.at(idStat) = nullptr;
                someDisabled = TRUE;
              }
            }
        }
        while(someDisabled);
      } // if(!fileNameNoNetPositions.empty())

      selectedNoNetStations = selectorNoNetStations->select(slr->times.front(), slr->times.back(), platforms);
      const UInt countStation = std::count(selectedNoNetStations.begin(), selectedNoNetStations.end(), TRUE);
      logInfo<<"  "<<countStation<<" stations contribute to the computation of net translation/rotation/scale"<<Log::endl;
      if(!countStation)
        throw(Exception("no stations contribute to the computation of net translation/rotation/scale"));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SlrParametrizationStaticPositions::initParameter(SlrNormalEquationInfo &normalEquationInfo)
{
  try
  {
    std::fill(index.begin(), index.end(), SlrParameterIndex());
    applyConstraint = FALSE;
    if(!isEnabled(normalEquationInfo, name))
      return;

    UInt countPara = 0;
    for(auto stat : slr->stations)
      if(stat->useable() && selectedStations.at(stat->idStat()) && normalEquationInfo.estimateStation.at(stat->idStat()))
      {
        std::vector<ParameterName> parameterNames({{stat->name(), "position.x"}, {stat->name(), "position.y"}, {stat->name(), "position.z"}});
        index.at(stat->idStat()) = normalEquationInfo.parameterNamesStation(stat->idStat(), parameterNames);
        countPara += parameterNames.size();
      }
    if(countPara)
      logInfo<<countPara%"%9i station static position parameters"s<<Log::endl;

    applyConstraint = isEnabled(normalEquationInfo, nameConstraint) && (sigmaNoNetRotation || sigmaNoNetTranslation) && countPara;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SlrParametrizationStaticPositions::aprioriParameter(const SlrNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const
{
  try
  {
    for(UInt idStat=0; idStat<index.size(); idStat++)
      if(index.at(idStat))
        copy(pos.at(idStat).vector(), x0.row(normalEquationInfo.index(index.at(idStat)), 3));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SlrParametrizationStaticPositions::designMatrix(const SlrNormalEquationInfo &/*normalEquationInfo*/, const SlrObservationEquation &eqn, SlrDesignMatrix &A) const
{
  try
  {
    if(index.at(eqn.station->idStat()))
    {
      MatrixSlice Design(A.column(index.at(eqn.station->idStat())));
      for(UInt idEpoch=0; idEpoch<eqn.timesTrans.size(); idEpoch++)
        matMult(1., eqn.A.slice(idEpoch, SlrObservationEquation::idxPosStat, 1, 3),
                slr->rotationCrf2Trf(eqn.timesTrans.at(idEpoch)).matrix().trans(),
                Design.row(idEpoch));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SlrParametrizationStaticPositions::constraints(const SlrNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const
{
  try
  {
    if(!applyConstraint)
      return;

    UInt noNetCount = 0;
    const UInt idxNNT = noNetCount; if(sigmaNoNetTranslation) noNetCount += 3;
    const UInt idxNNR = noNetCount; if(sigmaNoNetRotation)    noNetCount += 3;
    const UInt idxNNS = noNetCount; if(sigmaNoNetScale)       noNetCount += 1;

    UInt stationCount = 0;
    for(UInt idStat=0; idStat<slr->stations.size(); idStat++)
      if(selectedNoNetStations.at(idStat) && index.at(idStat))
        stationCount++;

    // observation equations
    Vector l(3*stationCount);
    Matrix A(l.rows(), noNetCount);
    UInt i = 0;
    for(UInt idStat=0; idStat<slr->stations.size(); idStat++)
      if(selectedNoNetStations.at(idStat) && index.at(idStat))
      {
        copy((noNetPos.at(idStat)-pos.at(idStat)).vector(), l.row(3*i, 3));
        if(sigmaNoNetTranslation)
          copy(identityMatrix(3), A.slice(3*i, idxNNT, 3, 3));
        if(sigmaNoNetRotation)
        {
          const Vector3d pos = noNetPos.at(idStat)/DEFAULT_R;
          A(3*i+0, idxNNR+1) =  pos.z(); A(3*i+0, idxNNR+2) = -pos.y();
          A(3*i+1, idxNNR+0) = -pos.z(); A(3*i+1, idxNNR+2) =  pos.x();
          A(3*i+2, idxNNR+0) =  pos.y(); A(3*i+2, idxNNR+1) = -pos.x();
        }
        if(sigmaNoNetScale)
        {
          const Vector3d posForScale = noNetPos.at(idStat)/DEFAULT_R;
          A(3*i+0, idxNNS+0) = posForScale.x();
          A(3*i+1, idxNNS+0) = posForScale.y();
          A(3*i+2, idxNNS+0) = posForScale.z();
        }
        i++;
      }

    // compute station weights for no-net constraint
    Vector x(A.columns()), sigma;
    if(rootMeanSquare(l) > 1e-4)
      x = Vce::robustLeastSquares(A, l, 3, huber, huberPower, 30, sigma);

    // compute noNetEstimator = (A'PA)^(-1) A'P
    for(UInt i=0; i<sigma.rows(); i++)
      A.row(3*i, 3) *= 1./sigma(i);
    const Vector tau = QR_decomposition(A);
    const Matrix W = A.row(0, noNetCount);
    generateQ(A, tau);
    triangularSolve(1., W, A.trans());
    for(UInt i=0; i<sigma.rows(); i++)
      A.row(3*i, 3) *= 1./sigma(i);
    noNetEstimator = A.trans();

    if(sigmaNoNetTranslation) logStatus<<"apply no-net translation to station positions, apriori ("<<1e3*x(idxNNT+0)%"%.1f, "s<<1e3*x(idxNNT+1)%"%.1f, "s<<1e3*x(idxNNT+2)%"%.1f) mm"s<<Log::endl;
    if(sigmaNoNetRotation)    logStatus<<"apply no-net rotation to station positions,    apriori ("<<1e3*x(idxNNR+0)%"%.1f, "s<<1e3*x(idxNNR+1)%"%.1f, "s<<1e3*x(idxNNR+2)%"%.1f) mm"s<<Log::endl;
    if(sigmaNoNetScale)       logStatus<<"apply no-net scale to receiver positions,      apriori ("<<1e3*x(idxNNS+0)%"%.1f) mm"s<<Log::endl;
    if(sigmaNoNetRotation || sigmaNoNetTranslation || sigmaNoNetScale)
    {
      logStatus<<"  no-net coordinate residuals rms = "<<1e3*rootMeanSquare(l-A*x)%"%.1f mm, "s<<Log::endl;
      UInt i = 0;
      if(sigma.size())
        for(UInt idStat=0; idStat<slr->stations.size(); idStat++)
          if(selectedNoNetStations.at(idStat) && index.at(idStat))
          {
            if(sigma(i) > std::pow(3./huber, huberPower))
              logWarning<<"  "<<slr->stations.at(idStat)->name()<<" outlier sigma = "<<sigma(i)%"%.2f"s<<Log::endl;
            i++;
          }
    }

    // weighted no-net constraints
    SlrDesignMatrix Design(normalEquationInfo, x);
    if(sigmaNoNetTranslation) Design.l.row(idxNNT, 3)  *= 1./sigmaNoNetTranslation;
    if(sigmaNoNetRotation)    Design.l.row(idxNNR, 3)  *= 1./sigmaNoNetRotation;
    if(sigmaNoNetScale)       Design.l.row(idxNNS, 1)  *= 1./sigmaNoNetScale;
    if(sigmaNoNetTranslation) A.trans().row(idxNNT, 3) *= 1./sigmaNoNetTranslation;
    if(sigmaNoNetRotation)    A.trans().row(idxNNR, 3) *= 1./sigmaNoNetRotation;
    if(sigmaNoNetScale)       A.trans().row(idxNNS, 1) *= 1./sigmaNoNetScale;
    i = 0;
    for(UInt idStat=0; idStat<slr->stations.size(); idStat++)
      if(selectedNoNetStations.at(idStat) && index.at(idStat))
        copy(A.row(3*i++,3).trans(), Design.column(index.at(idStat)));
    Design.accumulateNormals(normals, n, lPl, obsCount);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double SlrParametrizationStaticPositions::updateParameter(const SlrNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef /*Wz*/)
{
  try
  {
    if(!isEnabled(normalEquationInfo, name))
      return 0;

    // compute no-net translation and not-net rotation
    std::string infoNoNet;
    if(applyConstraint)
    {
      UInt noNetCount = 0;
      const UInt idxNNT = noNetCount; if(sigmaNoNetTranslation) noNetCount += 3;
      const UInt idxNNR = noNetCount; if(sigmaNoNetRotation)    noNetCount += 3;
      const UInt idxNNS = noNetCount; if(sigmaNoNetScale)       noNetCount += 1;
      Vector noNetPara(noNetCount);
      UInt i=0;
      for(UInt idStat=0; idStat<slr->stations.size(); idStat++)
        if(selectedNoNetStations.at(idStat) && index.at(idStat))
        {
          Vector dpos = (pos.at(idStat)-noNetPos.at(idStat)).vector();
          dpos += x.row(normalEquationInfo.index(index.at(idStat)), 3);
          matMult(1., noNetEstimator.column(3*i++, 3), dpos, noNetPara);
        }
      if(sigmaNoNetTranslation) infoNoNet += ", netTranslation ("+1e3*noNetPara(idxNNT+0)%"%.1f, "s+1e3*noNetPara(idxNNT+1)%"%.1f, "s+1e3*noNetPara(idxNNT+2)%"%.1f) mm"s;
      if(sigmaNoNetRotation)    infoNoNet += ", netRotation ("   +1e3*noNetPara(idxNNR+0)%"%.1f, "s+1e3*noNetPara(idxNNR+1)%"%.1f, "s+1e3*noNetPara(idxNNR+2)%"%.1f) mm"s;
      if(sigmaNoNetScale)       infoNoNet += ", netScale ("      +1e3*noNetPara(idxNNS+0)%"%.1f) mm"s;
    }

    // update positions
    Double maxChange = 0;
    Slr::InfoParameterChange info("mm");
    for(UInt idStat=0; idStat<index.size(); idStat++)
      if(index.at(idStat))
      {
        const Vector3d dpos(x.row(normalEquationInfo.index(index.at(idStat)), 3));
        pos.at(idStat) += dpos;
        slr->stations.at(idStat)->pos += dpos;
        if(info.update(1e3*dpos.r()))
          info.info = "position station ("+slr->stations.at(idStat)->name()+")"+infoNoNet;
      }
    info.print(1e-3, maxChange);
    return maxChange;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SlrParametrizationStaticPositions::writeResults(const SlrNormalEquationInfo &normalEquationInfo, const std::string &suffix) const
{
  try
  {
    if(!isEnabled(normalEquationInfo, name))
      return;

    if(!fileNameGrid.empty())
    {
      // collect positions
      Matrix data(slr->stations.size(), 6); // x, y, z , dx, dy, dz
      for(UInt idStat=0; idStat<slr->stations.size(); idStat++)
        if(index.at(idStat))
        {
          const Transform3d trf2local = slr->stations.at(idStat)->global2localFrame(slr->stations.at(idStat)->times.front());
          copy(pos.at(idStat).vector().trans(), data.slice(idStat, 0, 1, 3));
          copy(trf2local.transform(pos.at(idStat)-pos0.at(idStat)).vector().trans(), data.slice(idStat, 3, 1, 3)); // dpos in NEU (north, east, up) system
        }

      logStatus<<"write station grid file <"<<fileNameGrid.appendBaseName(suffix)<<">"<<Log::endl;
      std::vector<Vector3d>            point;
      std::vector<std::vector<Double>> value(3);
      for(UInt idStat=0; idStat<data.rows(); idStat++)
        if(index.at(idStat))
        {
          point.push_back(Vector3d(data.slice(idStat, 0, 1, 3)));
          value.at(0).push_back(data(idStat, 3)); // north
          value.at(1).push_back(data(idStat, 4)); // east
          value.at(2).push_back(data(idStat, 5)); // up
        }
      writeFileGriddedData(fileNameGrid.appendBaseName(suffix), GriddedData(Ellipsoid(), point, std::vector<Double>(point.size(), 1.), value));
    }

    if(!fileNamePosition.empty())
    {
      VariableList varList;
      varList.setVariable("station", "****");
      logStatus<<"write positions to files <"<<fileNamePosition(varList).appendBaseName(suffix)<<">"<<Log::endl;
      for(UInt idStat=0; idStat<index.size(); idStat++)
        if(index.at(idStat))
        {
          Vector3dEpoch epoch;
          epoch.time     = slr->stations.at(idStat)->times.at(static_cast<UInt>(std::floor(0.5*slr->stations.at(idStat)->times.size())));
          epoch.vector3d = pos.at(idStat);
          Vector3dArc arc;
          arc.push_back(epoch);
          varList.setVariable("station", slr->stations.at(idStat)->name());
          InstrumentFile::write(fileNamePosition(varList).appendBaseName(suffix), arc);
        }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
