/***********************************************/
/**
* @file observationMiscSstVariational.cpp
*
* @brief Satellite to satellite tracking (Variational equations).
*
* @author Torsten Mayer-Guerr
* @date 2015-06-02
*
*/
/***********************************************/

#include "base/import.h"
#include "files/fileInstrument.h"
#include "classes/parametrizationGravity/parametrizationGravity.h"
#include "classes/parametrizationAcceleration/parametrizationAcceleration.h"
#include "classes/parametrizationSatelliteTracking/parametrizationSatelliteTracking.h"
#include "misc/observation/variationalEquation.h"
#include "misc/observation/variationalEquationFromFile.h"
#include "observationMiscSstVariational.h"

/***********************************************/

ObservationMiscSstVariational::ObservationMiscSstVariational(Config &config)
{
  try
  {
    std::vector<FileName> fileNameSst;
    FileName              fileNamePod1, fileNamePod2;
    FileName              fileNameVariational1, fileNameVariational2;
    UInt                  integrationDegree;
    UInt                  interpolationDegree;
    sstType = 0;

    renameDeprecatedConfig(config, "representation", "parametrizationGravity",       date2time(2020, 6, 3));
    renameDeprecatedConfig(config, "parameter1",     "parametrizationAcceleration1", date2time(2020, 6, 3));
    renameDeprecatedConfig(config, "parameter2",     "parametrizationAcceleration2", date2time(2020, 6, 3));
    renameDeprecatedConfig(config, "parameterSst",   "parametrizationSst",           date2time(2020, 6, 3));

    if(readConfigSequence(config, "rightHandSide", Config::MUSTSET, "", "input for observation vectors"))
    {
      readConfig(config, "inputfileSatelliteTracking", fileNameSst,  Config::OPTIONAL, "", "ranging observations and corrections");
      readConfig(config, "inputfileOrbit1",            fileNamePod1, Config::OPTIONAL, "", "kinematic positions of satellite A as observations");
      readConfig(config, "inputfileOrbit2",            fileNamePod2, Config::OPTIONAL, "", "kinematic positions of satellite B as observations");
      endSequence(config);
    }
    std::string choice;
    if(readConfigChoice(config, "sstType", choice, Config::MUSTSET, "", ""))
    {
      if(readConfigChoiceElement(config, "range",     choice, "")) sstType = 0;
      if(readConfigChoiceElement(config, "rangeRate", choice, "")) sstType = 1;
      if(readConfigChoiceElement(config, "none",      choice, "")) sstType = 99;
      endChoice(config);
    }
    readConfig(config, "inputfileVariational1",        fileNameVariational1,   Config::MUSTSET,  "",    "approximate position and integrated state matrix");
    readConfig(config, "inputfileVariational2",        fileNameVariational2,   Config::MUSTSET,  "",    "approximate position and integrated state matrix");
    readConfig(config, "ephemerides",                  ephemerides,            Config::OPTIONAL, "jpl", "");
    readConfig(config, "parametrizationGravity",       parameterGravity,       Config::DEFAULT,  "",    "gravity field parametrization");
    readConfig(config, "parametrizationAcceleration1", parameterAcceleration1, Config::DEFAULT,  "",    "orbit1 force parameters");
    readConfig(config, "parametrizationAcceleration2", parameterAcceleration2, Config::DEFAULT,  "",    "orbit2 force parameters");
    readConfig(config, "parametrizationSst",           parameterSst,           Config::DEFAULT,  "",    "satellite tracking parameter");
    readConfig(config, "integrationDegree",            integrationDegree,      Config::DEFAULT,  "7",   "integration of forces by polynomial approximation of degree n");
    readConfig(config, "interpolationDegree",          interpolationDegree,    Config::DEFAULT,  "7",   "orbit interpolation by polynomial approximation of degree n");
    if(isCreateSchema(config)) return;

    // =======================

    // init
    // ----
    sstFile.resize(fileNameSst.size());
    for(UInt i=0; i<fileNameSst.size(); i++)
      sstFile.at(i) = InstrumentFile::newFile(fileNameSst.at(i));
    pod1File.open(fileNamePod1);
    pod2File.open(fileNamePod2);

    // test instrument files
    // ---------------------
    for(UInt i=1; i<sstFile.size(); i++)
      InstrumentFile::checkArcCount({*sstFile.at(i), *sstFile.at(0)});

    variationalEquation1.open(fileNameVariational1, parameterGravity, parameterAcceleration1, std::vector<Time>(), ephemerides, integrationDegree);
    variationalEquation2.open(fileNameVariational2, parameterGravity, parameterAcceleration2, std::vector<Time>(), ephemerides, integrationDegree);

    countArc = sstFile.at(0)->arcCount();
    computeVelocity = (sstType==1);
    polynomial.init(interpolationDegree);

    // =======================

    // count parameters
    // ----------------
    gravityCount = variationalEquation1.parameterCountGravity();
    state1Count  = variationalEquation1.parameterCount() - gravityCount;
    state2Count  = variationalEquation2.parameterCount() - gravityCount;

    countAParameter = 0;
    idxGravity  = countAParameter; countAParameter += gravityCount;
    idxState1   = countAParameter; countAParameter += state1Count;
    idxState2   = countAParameter; countAParameter += state2Count;
    idxSstPara  = countAParameter; countAParameter += parameterSst->parameterCount();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ObservationMiscSstVariational::setInterval(const Time &timeStart, const Time &timeEnd)
{
  try
  {
    parameterGravity->setInterval(timeStart, timeEnd);
    parameterAcceleration1->setInterval(timeStart, timeEnd);
    parameterAcceleration2->setInterval(timeStart, timeEnd);
    variationalEquation1.computeIndices();
    variationalEquation2.computeIndices();

    // count parameters
    // ----------------
    gravityCount = variationalEquation1.parameterCountGravity();
    state1Count  = variationalEquation1.parameterCount() - gravityCount;
    state2Count  = variationalEquation2.parameterCount() - gravityCount;

    countAParameter = 0;
    idxGravity  = countAParameter; countAParameter += gravityCount;
    idxState1   = countAParameter; countAParameter += state1Count;
    idxState2   = countAParameter; countAParameter += state2Count;
    idxSstPara  = countAParameter; countAParameter += parameterSst->parameterCount();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ObservationMiscSstVariational::parameterName(std::vector<ParameterName> &name) const
{
  try
  {
    parameterGravity->parameterName(name);

    std::vector<ParameterName> name1;
    variationalEquation1.parameterName(name1);
    const std::string satelliteName1 = variationalEquation1.satellite() ? variationalEquation1.satellite()->satelliteName : "satellite1";
    for(UInt i=gravityCount; i<name1.size(); i++)
    {
      name1.at(i).object = satelliteName1;
      name.push_back(name1.at(i));
    }

    std::vector<ParameterName> name2;
    variationalEquation2.parameterName(name2);
    const std::string satelliteName2 = variationalEquation2.satellite() ? variationalEquation2.satellite()->satelliteName : "satellite2";
    for(UInt i=gravityCount; i<name2.size(); i++)
    {
      name2.at(i).object = satelliteName2;
      name.push_back(name2.at(i));
    }

    if(parameterSst)
    {
      parameterSst->parameterName(name);
      for(UInt i=name.size()-parameterSst->parameterCount(); i<name.size(); i++)
        name.at(i).object = satelliteName1+"."+satelliteName2;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

ObservationMiscSstVariational::Arc ObservationMiscSstVariational::computeArc(UInt arcNo, CovarianceSstPtr covSst, CovariancePodPtr covPod1, CovariancePodPtr covPod2,
                                                                             const std::vector<Rotary3d> &rotSat1, const std::vector<Rotary3d> &rotSat2)
{
  try
  {
    // read SST observations
    // ---------------------
    std::vector<SatelliteTrackingArc> sst(sstFile.size());
    for(UInt k=0; k<sstFile.size(); k++)
      sst.at(k) = sstFile.at(k)->readArc(arcNo);
    for(UInt k=1; k<sst.size(); k++)
      ::Arc::checkSynchronized({sst.at(0), sst.at(k)});

    // read POD observations
    // ---------------------
    OrbitArc pod1 = pod1File.readArc(arcNo);
    OrbitArc pod2 = pod2File.readArc(arcNo);

    // =============================================

    const UInt epochCount = sst.at(0).size();
    parameterSst->setIntervalArc(sst.at(0).at(0).time, sst.at(0).back().time+medianSampling(sst.at(0).times()));

    // count observations and calculate index
    // --------------------------------------
    const UInt countSst  = (sst.size() && (sstType!=99)) ? sst.at(0).size() : 0;
    const UInt countPod1 = (pod1.size()<2) ? 0 : 3*pod1.size();
    const UInt countPod2 = (pod2.size()<2) ? 0 : 3*pod2.size();

    UInt obsCount = 0;
    const UInt idxSst   = obsCount; obsCount += countSst;
    const UInt idxPod1  = obsCount; obsCount += countPod1;
    const UInt idxPod2  = obsCount; obsCount += countPod2;

    // arc related parameters (in matrix B)
    // ------------------------------------
    UInt countBParameter = 0;
    const UInt idxSstParaArc = countBParameter; countBParameter += parameterSst->parameterCountArc();

    if(obsCount <= countBParameter)
      return Arc();

    Matrix l(obsCount, 1);
    Matrix A(obsCount, countAParameter);
    Matrix B(obsCount, countBParameter);

    // =============================================

    std::vector<Time> timesSst = sst.at(0).times();
    VariationalEquationFromFile::ObservationEquation eqn1 = variationalEquation1.integrateArc(timesSst.at(0), timesSst.back(), TRUE/*position*/, computeVelocity, rotSat1);
    VariationalEquationFromFile::ObservationEquation eqn2 = variationalEquation2.integrateArc(timesSst.at(0), timesSst.back(), TRUE/*position*/, computeVelocity, rotSat2);

    // =============================================

    // kinematic orbit observations
    // ----------------------------
    if(countPod1)
    {
      for(UInt k=0; k<countPod1/3; k++)
      {
        l(3*k+0+idxPod1,0) = pod1.at(k).position.x();
        l(3*k+1+idxPod1,0) = pod1.at(k).position.y();
        l(3*k+2+idxPod1,0) = pod1.at(k).position.z();
      }
      l.row(idxPod1, countPod1) -= polynomial.interpolate(pod1.times(), eqn1.times, eqn1.pos0, 3);

      Matrix D = polynomial.interpolate(pod1.times(), eqn1.times, eqn1.PosDesign, 3);
      if(gravityCount)
        copy(D.column(0, gravityCount),         A.slice(idxPod1, idxGravity, countPod1, gravityCount));
      copy(D.column(gravityCount, state1Count), A.slice(idxPod1, idxState1,  countPod1, state1Count));
    }

    if(countPod2)
    {
      for(UInt k=0; k<countPod2/3; k++)
      {
        l(3*k+0+idxPod2,0) = pod2.at(k).position.x();
        l(3*k+1+idxPod2,0) = pod2.at(k).position.y();
        l(3*k+2+idxPod2,0) = pod2.at(k).position.z();
      }
      l.row(idxPod2, countPod2) -= polynomial.interpolate(pod2.times(), eqn1.times, eqn2.pos0, 3);

      Matrix D = polynomial.interpolate(pod2.times(), eqn1.times, eqn2.PosDesign, 3);
      if(gravityCount)
        copy(D.column(0, gravityCount),         A.slice(idxPod2, idxGravity, countPod2, gravityCount));
      copy(D.column(gravityCount, state2Count), A.slice(idxPod2, idxState2,  countPod2, state2Count));
    }

    // =============================================

    if(countSst)
    {
      // inter satellite vectors
      // -----------------------
      std::vector<Vector3d> pos12(epochCount);
      std::vector<Vector3d> vel12(epochCount);
      std::vector<Vector3d> e12  (epochCount);
      Vector sst0(epochCount);

      Vector dpos0 = polynomial.interpolate(timesSst, eqn1.times, eqn2.pos0-eqn1.pos0, 3);
      for(UInt i=0; i<epochCount; i++)
      {
        pos12.at(i) = Vector3d(dpos0(3*i+0), dpos0(3*i+1), dpos0(3*i+2));
        e12.at(i)   = normalize(pos12.at(i));
      }

      if(computeVelocity)
      {
        Vector dvel0 = polynomial.interpolate(timesSst, eqn1.times, eqn2.vel0-eqn1.vel0, 3);
        for(UInt i=0; i<epochCount; i++)
          vel12.at(i) = Vector3d(dvel0(3*i+0), dvel0(3*i+1), dvel0(3*i+2));
      }

      // range
      if(sstType==0)
      {
        Matrix PosGravity;
        if(gravityCount) axpy(-1, eqn1.PosDesign.column(0, gravityCount), eqn2.PosDesign.column(0, gravityCount)); // gravity field difference
        if(gravityCount) PosGravity = polynomial.interpolate(timesSst, eqn2.times, eqn2.PosDesign.column(0, gravityCount));
        Matrix PosState1  = polynomial.interpolate(timesSst, eqn1.times, eqn1.PosDesign.column(gravityCount,state1Count), 3);
        Matrix PosState2  = polynomial.interpolate(timesSst, eqn2.times, eqn2.PosDesign.column(gravityCount,state2Count), 3);

        Double rho0  = pos12.at(0).r();
        for(UInt i=0; i<epochCount; i++)
        {
          Matrix e    = e12.at(i).vector().trans();
          Double rho  = pos12.at(i).r();
          sst0(i) = rho;

          for(UInt k=0; k<sst.size(); k++)
            l(idxSst+i,0) += (sst.at(k).at(i).range - sst.at(k).at(0).range);
          l(idxSst+i,0) -= (rho-rho0);

          // dRange/dPos12
          if(gravityCount)
            matMult( 1., e, PosGravity.row(3*i,3), A.slice(idxSst+i, idxGravity, 1, gravityCount));
          matMult(-1., e, PosState1.row(3*i,3),  A.slice(idxSst+i, idxState1,  1, state1Count));
          matMult( 1., e, PosState2.row(3*i,3),  A.slice(idxSst+i, idxState2,  1, state2Count));
        }
      } // if(sstType==0)

      // range rate
      if(sstType==1)
      {
        Matrix PosGravity, VelGravity;
        if(gravityCount) axpy(-1, eqn1.PosDesign.column(0, gravityCount), eqn2.PosDesign.column(0, gravityCount)); // gravity field difference
        if(gravityCount) axpy(-1, eqn1.VelDesign.column(0, gravityCount), eqn2.VelDesign.column(0, gravityCount)); // gravity field difference
        if(gravityCount) PosGravity = polynomial.interpolate(timesSst, eqn2.times, eqn2.PosDesign.column(0, gravityCount), 3);
        if(gravityCount) VelGravity = polynomial.interpolate(timesSst, eqn2.times, eqn2.VelDesign.column(0, gravityCount), 3);
        Matrix PosState1  = polynomial.interpolate(timesSst, eqn1.times, eqn1.PosDesign.column(gravityCount,state1Count), 3);
        Matrix PosState2  = polynomial.interpolate(timesSst, eqn2.times, eqn2.PosDesign.column(gravityCount,state2Count), 3);
        Matrix VelState1  = polynomial.interpolate(timesSst, eqn1.times, eqn1.VelDesign.column(gravityCount,state1Count), 3);
        Matrix VelState2  = polynomial.interpolate(timesSst, eqn2.times, eqn2.VelDesign.column(gravityCount,state2Count), 3);

        for(UInt i=0; i<epochCount; i++)
        {
          Double rho  = pos12.at(i).r();
          Double drho = inner(e12.at(i),vel12.at(i));

          for(UInt k=0; k<sst.size(); k++)
            l(idxSst+i,0) += sst.at(k).at(i).rangeRate;
          l(idxSst+i,0) -= drho;
          sst0(i) = drho;

          // dRate/dPos12
          Matrix e = ((1/rho)*vel12.at(i) - (drho/rho) * e12.at(i)).vector().trans();
          if(gravityCount)
            matMult( 1., e, PosGravity.row(3*i,3), A.slice(idxSst+i, idxGravity, 1, gravityCount));
          matMult(-1., e, PosState1.row(3*i,3),  A.slice(idxSst+i, idxState1,  1, state1Count));
          matMult( 1., e, PosState2.row(3*i,3),  A.slice(idxSst+i, idxState2,  1, state2Count));

          // dRate/dVel12 = e12
          e = e12.at(i).vector().trans();
          if(gravityCount)
            matMult( 1., e, VelGravity.row(3*i,3), A.slice(idxSst+i, idxGravity, 1, gravityCount));
          matMult(-1., e, VelState1.row(3*i,3),  A.slice(idxSst+i, idxState1,  1, state1Count));
          matMult( 1., e, VelState2.row(3*i,3),  A.slice(idxSst+i, idxState2,  1, state2Count));
        }
      } // if(sstType==1)

      if(parameterSst)
      {
        parameterSst->compute(sstType, timesSst, sst0, eqn1.pos0, eqn2.pos0, eqn1.vel0, eqn2.vel0, eqn1.rotSat, eqn2.rotSat,
                              A.slice(idxSst, idxSstPara,    countSst, parameterSst->parameterCount()),
                              B.slice(idxSst, idxSstParaArc, countSst, parameterSst->parameterCountArc()));
      }
    } // if(countSst)

    // =============================================

    // decorrelation
    // -------------
    if(covSst  && countSst)  covSst->decorrelate (arcNo, timesSst, {l.row(idxSst,  countSst),  A.row(idxSst,  countSst),  B.row(idxSst,  countSst)});
    if(covPod1 && countPod1) covPod1->decorrelate(arcNo, pod1,     {l.row(idxPod1, countPod1), A.row(idxPod1, countPod1), B.row(idxPod1, countPod1)});
    if(covPod2 && countPod2) covPod2->decorrelate(arcNo, pod2,     {l.row(idxPod2, countPod2), A.row(idxPod2, countPod2), B.row(idxPod2, countPod2)});

   // =============================================

    Arc observationArc;
    observationArc.l = l;
    observationArc.A = A;
    observationArc.B = B;
    observationArc.timesSst  = timesSst;
    observationArc.timesPod1 = pod1.times();
    observationArc.timesPod2 = pod2.times();
    observationArc.pod1      = pod1;
    observationArc.pod2      = pod2;
    observationArc.rotSat1   = eqn1.rotSat;
    observationArc.rotSat2   = eqn2.rotSat;
    observationArc.pos1      = eqn1.pos0;
    observationArc.pos2      = eqn2.pos0;
    return observationArc;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/
