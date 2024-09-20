/***********************************************/
/**
* @file gnssParametrizationStaticPositions.cpp
*
* @brief Position estimation with no-net constraints.
* @see GnssParametrization
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @author HongzhanZhao
* @date 2021-01-23
*
*/
/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "files/fileInstrument.h"
#include "files/fileGriddedData.h"
#include "classes/platformSelector/platformSelector.h"
#include "misc/varianceComponentEstimation.h"
#include "gnss/gnssParametrization/gnssParametrization.h"
#include "gnss/gnssParametrization/gnssParametrizationStaticPositions.h"

/***********************************************/

GnssParametrizationStaticPositions::GnssParametrizationStaticPositions(Config &config)
{
  try
  {
    sigmaNoNetTranslation = sigmaNoNetRotation = sigmaNoNetScale = 0;
    huber = 1e99; huberPower = 0;

    readConfig(config, "name",                      name,                   Config::OPTIONAL, "parameter.staticPositions", "used for parameter selection");
    readConfig(config, "selectReceivers",           selectReceivers,        Config::MUSTSET,  "", "");
    readConfig(config, "outputfileGriddedPosition", fileNameGrid,           Config::OPTIONAL, "output/gridPosition_{loopTime:%D}.dat",              "delta north east up for all stations");
    readConfig(config, "outputfilePosition",        fileNamePosition,       Config::OPTIONAL, "output/stationPosition_{loopTime:%D}.{station}.dat", "variable {station} available, full estimated coordinates (in TRF)");
    readConfig(config, "nameConstraint",            nameConstraint,         Config::OPTIONAL, "constraint.staticPositions", "used for parameter selection");
    readConfig(config, "selectNoNetReceivers",      selectNoNetReceivers,   Config::DEFAULT,  R"(["all"])", "");
    readConfig(config, "inputfileNoNetPositions",   fileNameNoNetPositions, Config::OPTIONAL, "{groopsDataDir}/gnss/receiverStation/position/igs/igs20/stationPosition.{station}.dat", "variable {station} available, precise coordinates used for no-net constraints (in TRF)");
    readConfig(config, "noNetTranslationSigma",     sigmaNoNetTranslation,  Config::OPTIONAL, "0.01",  "(0 = unconstrained) sigma [m] for no-net translation constraint on station coordinates");
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

void GnssParametrizationStaticPositions::init(Gnss *gnss, Parallel::CommunicatorPtr comm)
{
  try
  {
    this->gnss  = gnss;
    index.resize(gnss->receivers.size());
    selectedReceivers = gnss->selectReceivers(selectReceivers);

    // apriori positions
    pos0.resize(gnss->receivers.size());
    for(auto recv : gnss->receivers)
      if(selectedReceivers.at(recv->idRecv()) && recv->useable())
      {
        if(!recv->isEarthFixed())
          throw(Exception(recv->name()+" must be given in an Earth fixed frame"));
        Vector p(3);
        if(recv->isMyRank())
        {
          UInt count = 0;
          for(UInt idEpoch=0; idEpoch<gnss->times.size(); idEpoch++)
            if(recv->useable(idEpoch))
            {
              p += recv->pos.at(idEpoch).vector();
              count++;
            }
          p *= 1./count;
        }
        Parallel::reduceSum(p, 0, comm);
        Parallel::broadCast(p, 0, comm);
        pos0.at(recv->idRecv()) = Vector3d(p);
      }

    // no net-positions
    // ----------------
    if(sigmaNoNetRotation || sigmaNoNetTranslation || sigmaNoNetScale)
    {
      std::vector<const Platform*> platforms(gnss->receivers.size(), nullptr);
      for(auto recv : gnss->receivers)
        if(selectedReceivers.at(recv->idRecv()) && recv->useable())
          platforms.at(recv->idRecv()) = &recv->platform;

      // no net positions
      noNetPos = pos0;
      if(!fileNameNoNetPositions.empty())
        for(UInt idRecv=0; idRecv<gnss->receivers.size(); idRecv++)
          if(platforms.at(idRecv))
          {
            try
            {
              VariableList fileNameVariableList;
              fileNameVariableList.setVariable("station", gnss->receivers.at(idRecv)->name());
              Vector3dArc arc = InstrumentFile::read(fileNameNoNetPositions(fileNameVariableList));
              auto iter = (arc.size() == 1) ? arc.begin() : std::find_if(arc.begin(), arc.end(), [&](const Epoch &e){return e.time.isInInterval(gnss->times.front(), gnss->times.back());});
              if(iter == arc.end())
                throw(Exception("no position found"));
              noNetPos.at(idRecv) = iter->vector3d;
            }
            catch(std::exception &/*e*/)
            {
              platforms.at(idRecv) = nullptr;
            }
          }

      selectedNoNetReceivers = selectNoNetReceivers->select(gnss->times.front(), gnss->times.back(), platforms);
      const UInt countStation = std::count(selectedNoNetReceivers.begin(), selectedNoNetReceivers.end(), TRUE);
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

void GnssParametrizationStaticPositions::initParameter(GnssNormalEquationInfo &normalEquationInfo)
{
  try
  {
    std::fill(index.begin(), index.end(), GnssParameterIndex());
    applyConstraint = FALSE;
    if(!isEnabled(normalEquationInfo, name))
      return;

    UInt countPara = 0;
    for(auto recv : gnss->receivers)
      if(recv->useable() && selectedReceivers.at(recv->idRecv()) && normalEquationInfo.estimateReceiver.at(recv->idRecv()))
      {
        std::vector<ParameterName> parameterNames({{recv->name(), "position.x"}, {recv->name(), "position.y"}, {recv->name(), "position.z"}});
        index.at(recv->idRecv()) = normalEquationInfo.parameterNamesReceiver(recv->idRecv(), parameterNames);
        countPara += parameterNames.size();
      }
    if(countPara)
      logInfo<<countPara%"%9i receiver static position parameters"s<<Log::endl;

    applyConstraint = isEnabled(normalEquationInfo, nameConstraint) && !normalEquationInfo.isEachReceiverSeparately
                    && (sigmaNoNetRotation || sigmaNoNetTranslation || sigmaNoNetScale) && countPara;

    // synchronize positions
    pos.resize(gnss->receivers.size());
    for(auto recv : gnss->receivers)
      if(index.at(recv->idRecv()))
      {
        Vector p(3);
        if(recv->isMyRank())
        {
          UInt count = 0;
          for(UInt idEpoch=0; idEpoch<gnss->times.size(); idEpoch++)
            if(recv->useable(idEpoch))
            {
              p += recv->pos.at(idEpoch).vector();
              count++;
            }
          p *= 1./count;
        }
        Parallel::reduceSum(p, 0, normalEquationInfo.comm);
        Parallel::broadCast(p, 0, normalEquationInfo.comm);
        pos.at(recv->idRecv()) = Vector3d(p);
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationStaticPositions::aprioriParameter(const GnssNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const
{
  try
  {
    if(Parallel::isMaster(normalEquationInfo.comm))
      for(UInt idRecv=0; idRecv<index.size(); idRecv++)
        if(index.at(idRecv))
          copy(pos.at(idRecv).vector(), x0.row(normalEquationInfo.index(index.at(idRecv)), 3));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationStaticPositions::designMatrix(const GnssNormalEquationInfo &/*normalEquationInfo*/, const GnssObservationEquation &eqn, GnssDesignMatrix &A) const
{
  try
  {
    if(index.at(eqn.receiver->idRecv()))
      matMult(1., eqn.A.column(GnssObservationEquation::idxPosRecv, 3),
              gnss->rotationCrf2Trf(eqn.timeRecv).matrix().trans(),
              A.column(index.at(eqn.receiver->idRecv())));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationStaticPositions::constraints(const GnssNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const
{
  try
  {
    if(!Parallel::isMaster(normalEquationInfo.comm) || !applyConstraint)
      return;

    UInt noNetCount = 0;
    const UInt idxNNT = noNetCount; if(sigmaNoNetTranslation) noNetCount += 3;
    const UInt idxNNR = noNetCount; if(sigmaNoNetRotation)    noNetCount += 3;
    const UInt idxNNS = noNetCount; if(sigmaNoNetScale)       noNetCount += 1;

    UInt stationCount = 0;
    for(UInt idRecv=0; idRecv<gnss->receivers.size(); idRecv++)
      if(selectedNoNetReceivers.at(idRecv) && index.at(idRecv))
        stationCount++;

    // observation equations
    Vector l(3*stationCount);
    Matrix A(l.rows(), noNetCount);
    UInt i = 0;
    for(UInt idRecv=0; idRecv<gnss->receivers.size(); idRecv++)
      if(selectedNoNetReceivers.at(idRecv) && index.at(idRecv))
      {
        copy((noNetPos.at(idRecv)-pos.at(idRecv)).vector(), l.row(3*i, 3));
        if(sigmaNoNetTranslation)
          copy(identityMatrix(3), A.slice(3*i, idxNNT, 3, 3));
        if(sigmaNoNetRotation)
        {
          const Vector3d pos = noNetPos.at(idRecv)/DEFAULT_R;
          A(3*i+0, idxNNR+1) =  pos.z(); A(3*i+0, idxNNR+2) = -pos.y();
          A(3*i+1, idxNNR+0) = -pos.z(); A(3*i+1, idxNNR+2) =  pos.x();
          A(3*i+2, idxNNR+0) =  pos.y(); A(3*i+2, idxNNR+1) = -pos.x();
        }
        if(sigmaNoNetScale)
        {
          const Vector3d posForScale = noNetPos.at(idRecv)/DEFAULT_R;
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

    if(sigmaNoNetTranslation) logStatus<<"apply no-net translation to receiver positions, apriori ("<<1e3*x(idxNNT+0)%"%.1f, "s<<1e3*x(idxNNT+1)%"%.1f, "s<<1e3*x(idxNNT+2)%"%.1f) mm"s<<Log::endl;
    if(sigmaNoNetRotation)    logStatus<<"apply no-net rotation to receiver positions,    apriori ("<<1e3*x(idxNNR+0)%"%.1f, "s<<1e3*x(idxNNR+1)%"%.1f, "s<<1e3*x(idxNNR+2)%"%.1f) mm"s<<Log::endl;
    if(sigmaNoNetScale)       logStatus<<"apply no-net scale to receiver positions,       apriori ("<<1e3*x(idxNNS+0)%"%.1f) mm"s<<Log::endl;
    if(sigmaNoNetRotation || sigmaNoNetTranslation || sigmaNoNetScale)
    {
      logStatus<<"  no-net coordinate residuals rms = "<<1e3*rootMeanSquare(l-A*x)%"%.1f mm, "s<<Log::endl;
      UInt i = 0;
      if(sigma.size())
        for(UInt idRecv=0; idRecv<gnss->receivers.size(); idRecv++)
          if(selectedNoNetReceivers.at(idRecv) && index.at(idRecv))
          {
            if(sigma(i) > std::pow(3./huber, huberPower))
              logWarning<<"  "<<gnss->receivers.at(idRecv)->name()<<" outlier sigma = "<<sigma(i)%"%.2f"s<<Log::endl;
            i++;
          }
    }

    // weighted no-net constraints
    GnssDesignMatrix Design(normalEquationInfo, x);
    if(sigmaNoNetTranslation) Design.l.row(idxNNT, 3)  *= 1./sigmaNoNetTranslation;
    if(sigmaNoNetRotation)    Design.l.row(idxNNR, 3)  *= 1./sigmaNoNetRotation;
    if(sigmaNoNetScale)       Design.l.row(idxNNS, 1)  *= 1./sigmaNoNetScale;
    if(sigmaNoNetTranslation) A.trans().row(idxNNT, 3) *= 1./sigmaNoNetTranslation;
    if(sigmaNoNetRotation)    A.trans().row(idxNNR, 3) *= 1./sigmaNoNetRotation;
    if(sigmaNoNetScale)       A.trans().row(idxNNS, 1) *= 1./sigmaNoNetScale;
    i = 0;
    for(UInt idRecv=0; idRecv<gnss->receivers.size(); idRecv++)
      if(selectedNoNetReceivers.at(idRecv) && index.at(idRecv))
        copy(A.row(3*i++,3).trans(), Design.column(index.at(idRecv)));
    Design.accumulateNormals(normals, n, lPl, obsCount);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double GnssParametrizationStaticPositions::updateParameter(const GnssNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef /*Wz*/)
{
  try
  {
    if(!isEnabled(normalEquationInfo, name))
      return 0;

    // compute no-net translation and not-net rotation
    std::string infoNoNet;
    if(applyConstraint)
    {
      Parallel::broadCast(noNetEstimator, 0, normalEquationInfo.comm);
      UInt noNetCount = 0;
      const UInt idxNNT = noNetCount; if(sigmaNoNetTranslation) noNetCount += 3;
      const UInt idxNNR = noNetCount; if(sigmaNoNetRotation)    noNetCount += 3;
      const UInt idxNNS = noNetCount; if(sigmaNoNetScale)       noNetCount += 1;
      Vector noNetPara(noNetCount);
      UInt i=0;
      for(UInt idRecv=0; idRecv<gnss->receivers.size(); idRecv++)
        if(selectedNoNetReceivers.at(idRecv) && index.at(idRecv))
        {
          Vector dpos = (pos.at(idRecv)-noNetPos.at(idRecv)).vector();
          dpos += x.row(normalEquationInfo.index(index.at(idRecv)), 3);
          matMult(1., noNetEstimator.column(3*i++, 3), dpos, noNetPara);
        }
      if(sigmaNoNetTranslation) infoNoNet += ", netTranslation ("+1e3*noNetPara(idxNNT+0)%"%.1f, "s+1e3*noNetPara(idxNNT+1)%"%.1f, "s+1e3*noNetPara(idxNNT+2)%"%.1f) mm"s;
      if(sigmaNoNetRotation)    infoNoNet += ", netRotation ("   +1e3*noNetPara(idxNNR+0)%"%.1f, "s+1e3*noNetPara(idxNNR+1)%"%.1f, "s+1e3*noNetPara(idxNNR+2)%"%.1f) mm"s;
      if(sigmaNoNetScale)       infoNoNet += ", netScale ("      +1e3*noNetPara(idxNNS+0)%"%.1f) mm"s;
    }

    // update positions
    Double maxChange = 0;
    Gnss::InfoParameterChange info("mm");
    for(UInt idRecv=0; idRecv<index.size(); idRecv++)
      if(index.at(idRecv))
      {
        const Vector3d dpos(x.row(normalEquationInfo.index(index.at(idRecv)), 3));
        pos.at(idRecv) += dpos;
        if(gnss->receivers.at(idRecv)->isMyRank())
          for(UInt idEpoch=0; idEpoch<gnss->receivers.at(idRecv)->times.size(); idEpoch++)
            gnss->receivers.at(idRecv)->pos.at(idEpoch) += dpos;
        if(info.update(1e3*dpos.r()))
          info.info = "position receiver ("+gnss->receivers.at(idRecv)->name()+")"+infoNoNet;
      }
    info.synchronizeAndPrint(normalEquationInfo.comm, 1e-3, maxChange);
    return maxChange;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationStaticPositions::writeResults(const GnssNormalEquationInfo &normalEquationInfo, const std::string &suffix) const
{
  try
  {
    if(!isEnabled(normalEquationInfo, name))
      return;

    if(!fileNameGrid.empty())
    {
      // collect positions
      Matrix data(gnss->receivers.size(), 6); // x, y, z , dx, dy, dz
      for(UInt idRecv=0; idRecv<gnss->receivers.size(); idRecv++)
        if(index.at(idRecv) && gnss->receivers.at(idRecv)->isMyRank())
        {
          // find first valid Epoch
          Transform3d trf2local;
          for(UInt idEpoch : normalEquationInfo.idEpochs)
            if(gnss->receivers.at(idRecv)->useable(idEpoch))
            {
              trf2local = gnss->receivers.at(idRecv)->global2localFrame(idEpoch);
              break;
            }
          copy(pos.at(idRecv).vector().trans(), data.slice(idRecv, 0, 1, 3));
          copy(trf2local.transform(pos.at(idRecv)-pos0.at(idRecv)).vector().trans(), data.slice(idRecv, 3, 1, 3)); // dpos in NEU (north, east, up) system
        }
      Parallel::reduceSum(data, 0, normalEquationInfo.comm);

      if(Parallel::isMaster(normalEquationInfo.comm))
      {
        logStatus<<"write receiver grid file <"<<fileNameGrid.appendBaseName(suffix)<<">"<<Log::endl;
        std::vector<Vector3d>            point;
        std::vector<std::vector<Double>> value(3);
        for(UInt idRecv=0; idRecv<data.rows(); idRecv++)
          if(index.at(idRecv))
          {
            point.push_back(Vector3d(data.slice(idRecv, 0, 1, 3)));
            value.at(0).push_back(data(idRecv, 3)); // north
            value.at(1).push_back(data(idRecv, 4)); // east
            value.at(2).push_back(data(idRecv, 5)); // up
          }
        writeFileGriddedData(fileNameGrid.appendBaseName(suffix), GriddedData(Ellipsoid(), point, std::vector<Double>(point.size(), 1.), value));
      }
    }

    if(!fileNamePosition.empty())
    {
      VariableList fileNameVariableList;
      fileNameVariableList.setVariable("station", "****");
      logStatus<<"write positions to files <"<<fileNamePosition(fileNameVariableList).appendBaseName(suffix)<<">"<<Log::endl;
      for(UInt idRecv=0; idRecv<index.size(); idRecv++)
        if(index.at(idRecv) && gnss->receivers.at(idRecv)->isMyRank())
        {
          Vector3dEpoch epoch;
          epoch.time     = gnss->times.at(static_cast<UInt>(std::floor(0.5*gnss->times.size())));
          epoch.vector3d = pos.at(idRecv);
          Vector3dArc arc;
          arc.push_back(epoch);
          fileNameVariableList.setVariable("station", gnss->receivers.at(idRecv)->name());
          InstrumentFile::write(fileNamePosition(fileNameVariableList).appendBaseName(suffix), arc);
        }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
