/***********************************************/
/**
* @file observationMiscSstIntegral.cpp
*
* @brief Satellite to satellite tracking (Short Arc Integral).
*
* @author Torsten Mayer-Guerr
* @date 2015-06-02
*
*/
/***********************************************/

#include "base/import.h"
#include "files/fileInstrument.h"
#include "files/fileSatelliteModel.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/gravityfield/gravityfield.h"
#include "classes/parametrizationGravity/parametrizationGravity.h"
#include "classes/parametrizationAcceleration/parametrizationAcceleration.h"
#include "classes/parametrizationSatelliteTracking/parametrizationSatelliteTracking.h"
#include "misc/observation/integralEquation.h"
#include "misc/observation/covarianceSst.h"
#include "misc/observation/covariancePod.h"
#include "misc/observation/observationMisc.h"
#include "observationMiscSstIntegral.h"

/***********************************************/

ObservationMiscSstIntegral::ObservationMiscSstIntegral(Config &config)
{
  try
  {
    sstType = 0;
    std::string choice;
    FileName    fileNameSatellite1, fileNameSatellite2;
    FileName    orbit1Name, orbit2Name;
    FileName    starCamera1Name, starCamera2Name;

    renameDeprecatedConfig(config, "satelliteModel1", "inputfileSatelliteModel1",     date2time(2020, 8, 19));
    renameDeprecatedConfig(config, "satelliteModel2", "inputfileSatelliteModel2",     date2time(2020, 8, 19));
    renameDeprecatedConfig(config, "representation",  "parametrizationGravity",       date2time(2020, 6,  3));
    renameDeprecatedConfig(config, "parameter1",      "parametrizationAcceleration1", date2time(2020, 6,  3));
    renameDeprecatedConfig(config, "parameter2",      "parametrizationAcceleration2", date2time(2020, 6,  3));
    renameDeprecatedConfig(config, "parameterSst",    "parametrizationSst",           date2time(2020, 6,  3));

    readConfig(config, "inputfileSatelliteModel1", fileNameSatellite1,  Config::OPTIONAL, "{groopsDataDir}/satelliteModel/", "satellite macro model");
    readConfig(config, "inputfileSatelliteModel2", fileNameSatellite2,  Config::OPTIONAL, "{groopsDataDir}/satelliteModel/", "satellite macro model");
    readConfig(config, "rightHandSide",            rhs,                 Config::MUSTSET,  "", "input for the reduced observation vector");
    if(readConfigChoice(config, "sstType", choice, Config::MUSTSET, "", ""))
    {
      if(readConfigChoiceElement(config, "range",             choice, "")) sstType = 0;
      if(readConfigChoiceElement(config, "rangeRate",         choice, "")) sstType = 1;
      if(readConfigChoiceElement(config, "rangeAcceleration", choice, "")) sstType = 2;
      if(readConfigChoiceElement(config, "none",              choice, "")) sstType = 99;
      endChoice(config);
    }
    readConfig(config, "inputfileOrbit1",              orbit1Name,             Config::MUSTSET,  "",    "used to evaluate the observation equations, not used as observations");
    readConfig(config, "inputfileOrbit2",              orbit2Name,             Config::MUSTSET,  "",    "used to evaluate the observation equations, not used as observations");
    readConfig(config, "inputfileStarCamera1",         starCamera1Name,        Config::MUSTSET,  "",    "");
    readConfig(config, "inputfileStarCamera2",         starCamera2Name,        Config::MUSTSET,  "",    "");
    readConfig(config, "earthRotation",                earthRotation,          Config::MUSTSET,  "",    "");
    readConfig(config, "ephemerides",                  ephemerides,            Config::OPTIONAL, "jpl", "");
    readConfig(config, "gradientfield",                gradientfield,          Config::DEFAULT,  "",    "low order field to estimate the change of the gravity by position adjustement");
    readConfig(config, "parametrizationGravity",       parameterGravity,       Config::DEFAULT,  "",    "gravity field parametrization");
    readConfig(config, "parametrizationAcceleration1", parameterAcceleration1, Config::DEFAULT,  "",    "orbit1 force parameters");
    readConfig(config, "parametrizationAcceleration2", parameterAcceleration2, Config::DEFAULT,  "",    "orbit2 force parameters");
    readConfig(config, "parametrizationSst",           parameterSst,           Config::DEFAULT,  "",    "satellite tracking parameter");
    readConfig(config, "keepSatelliteStates",          keepSatelliteStates,    Config::DEFAULT,  "0",   "set boundary values of each arc global");    // satellite tracking type
    readConfig(config, "integrationDegree",            integrationDegree,      Config::DEFAULT,  "7",   "integration of forces by polynomial approximation of degree n");
    readConfig(config, "interpolationDegree",          interpolationDegree,    Config::DEFAULT,  "7",   "orbit interpolation by polynomial approximation of degree n");
    if(isCreateSchema(config)) return;

    // =======================

    if(!fileNameSatellite1.empty()) readFileSatelliteModel(fileNameSatellite1, satellite1);
    if(!fileNameSatellite2.empty()) readFileSatelliteModel(fileNameSatellite2, satellite2);

    // =======================

    // test instrument files
    // ---------------------
    orbit1File.open(orbit1Name);
    orbit2File.open(orbit2Name);
    starCamera1File.open(starCamera1Name);
    starCamera2File.open(starCamera2Name);

    std::vector<std::reference_wrapper<const InstrumentFile>> fileList{orbit1File, orbit2File, starCamera1File, starCamera2File};
    for(UInt j=0; j<rhs.size(); j++)
    {
      for(UInt i=0; i<rhs.at(j)->sstFile.size(); i++)
        fileList.push_back(*rhs.at(j)->sstFile.at(i));
      fileList.push_back(*rhs.at(j)->accelerometer1File);
      fileList.push_back(*rhs.at(j)->accelerometer2File);
    }
    InstrumentFile::checkArcCount(fileList);

    for(UInt j=1; j<rhs.size(); j++)
      InstrumentFile::checkArcCount({*rhs.at(j)->orbit1File, *rhs.at(0)->orbit1File});
    for(UInt j=1; j<rhs.size(); j++)
      InstrumentFile::checkArcCount({*rhs.at(j)->orbit2File, *rhs.at(0)->orbit2File});

    // init
    // ----
    countArc = orbit1File.arcCount();
    integralEquation.init(integrationDegree, interpolationDegree);

    computeRange        = (sstType==0);
    computeVelocity     = (sstType==1);
    computeAcceleration = (sstType==2);

    // =======================

    // design matrix A
    // ---------------
    countAParameter = 0;
    gravityCount = parameterGravity->parameterCount();
    idxGravity   = countAParameter; countAParameter += gravityCount;
    idxSat1      = countAParameter; countAParameter += parameterAcceleration1->parameterCount();
    idxSat2      = countAParameter; countAParameter += parameterAcceleration2->parameterCount();
    idxSstPara   = countAParameter; countAParameter += parameterSst->parameterCount();
    if(keepSatelliteStates)
    {
      idxBound1 = countAParameter; countAParameter += 6*countArc; // 2 boundary pos. (x,y,z).
      idxBound2 = countAParameter; countAParameter += 6*countArc; // 2 boundary pos. (x,y,z).
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool ObservationMiscSstIntegral::setInterval(const Time &timeStart, const Time &timeEnd)
{
  try
  {
    Bool change = FALSE;
    change = parameterGravity->setInterval(timeStart, timeEnd)       || change;
    change = parameterAcceleration1->setInterval(timeStart, timeEnd) || change;
    change = parameterAcceleration2->setInterval(timeStart, timeEnd) || change;
    change = parameterSst->setInterval(timeStart, timeEnd)           || change;
    if(!change)
      return FALSE;

    // count parameters
    // ----------------
    countAParameter = 0;
    gravityCount = parameterGravity->parameterCount();
    idxGravity   = countAParameter; countAParameter += gravityCount;
    idxSat1      = countAParameter; countAParameter += parameterAcceleration1->parameterCount();
    idxSat2      = countAParameter; countAParameter += parameterAcceleration2->parameterCount();
    idxSstPara   = countAParameter; countAParameter += parameterSst->parameterCount();
    if(keepSatelliteStates)
    {
      idxBound1 = countAParameter; countAParameter += 6*countArc; // 2 boundary pos. (x,y,z).
      idxBound2 = countAParameter; countAParameter += 6*countArc; // 2 boundary pos. (x,y,z).
    }
    return change;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ObservationMiscSstIntegral::parameterName(std::vector<ParameterName> &name) const
{
  try
  {
    parameterGravity->parameterName(name);

    const std::string satelliteName1 = satellite1 ? satellite1->satelliteName : "satellite1";
    const std::string satelliteName2 = satellite2 ? satellite2->satelliteName : "satellite2";

    parameterAcceleration1->parameterName(name);
    for(UInt i=name.size()-parameterAcceleration1->parameterCount(); i<name.size(); i++)
      name.at(i).object = satelliteName1;

    parameterAcceleration2->parameterName(name);
    for(UInt i=name.size()-parameterAcceleration2->parameterCount(); i<name.size(); i++)
      name.at(i).object = satelliteName2;

    if(parameterSst)
    {
      parameterSst->parameterName(name);
      for(UInt i=name.size()-parameterSst->parameterCount(); i<name.size(); i++)
        name.at(i).object = satelliteName1+"."+satelliteName2;
    }

    if(keepSatelliteStates)
    {
      for(UInt arcNo=0; arcNo<arcCount(); arcNo++)
      {
        const std::string str = "arc"+arcNo%"%i"s+"position.";
        name.push_back(ParameterName(satelliteName1, str+"start.x"));
        name.push_back(ParameterName(satelliteName1, str+"start.y"));
        name.push_back(ParameterName(satelliteName1, str+"start.z"));
        name.push_back(ParameterName(satelliteName1, str+"end.x"));
        name.push_back(ParameterName(satelliteName1, str+"end.y"));
        name.push_back(ParameterName(satelliteName1, str+"end.z"));
      }
      for(UInt arcNo=0; arcNo<arcCount(); arcNo++)
      {
        const std::string str = "arc"+arcNo%"%i"s+"position.";
        name.push_back(ParameterName(satelliteName2, str+"start.x"));
        name.push_back(ParameterName(satelliteName2, str+"start.y"));
        name.push_back(ParameterName(satelliteName2, str+"start.z"));
        name.push_back(ParameterName(satelliteName2, str+"end.x"));
        name.push_back(ParameterName(satelliteName2, str+"end.y"));
        name.push_back(ParameterName(satelliteName2, str+"end.z"));
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

ObservationMiscSst::Arc ObservationMiscSstIntegral::computeArc(UInt arcNo, CovarianceSstPtr covSst, CovariancePodPtr covPod1, CovariancePodPtr covPod2)
{
  try
  {
    OrbitArc      orbit1      = orbit1File.readArc(arcNo);
    OrbitArc      orbit2      = orbit2File.readArc(arcNo);
    StarCameraArc starCamera1 = starCamera1File.readArc(arcNo);
    StarCameraArc starCamera2 = starCamera2File.readArc(arcNo);
    ::Arc::checkSynchronized({orbit1, orbit2, starCamera1, starCamera2});
    const UInt    rhsCount    = rhs.size();
    const UInt    epochCount  = orbit1.size();

    parameterAcceleration1->setIntervalArc(orbit1.at(0).time, orbit1.back().time+medianSampling(orbit1.times()));
    parameterAcceleration2->setIntervalArc(orbit2.at(0).time, orbit2.back().time+medianSampling(orbit2.times()));
    parameterSst->setIntervalArc(orbit1.at(0).time, orbit1.back().time+medianSampling(orbit1.times()));

    // =============================================

    // count arc related parameters (in matrix B)
    // ------------------------------------
    UInt countBParameter = 0;
    if(!keepSatelliteStates)
    {
      idxBound1 = countBParameter; countBParameter += 6; // 2 Randpos. je (x,y,z).
      idxBound2 = countBParameter; countBParameter += 6; // 2 Randpos. je (x,y,z).
    }
    const UInt idxSat1Arc    = countBParameter; countBParameter += parameterAcceleration1->parameterCountArc();
    const UInt idxSat2Arc    = countBParameter; countBParameter += parameterAcceleration2->parameterCountArc();
    const UInt idxSstParaArc = countBParameter; countBParameter += parameterSst->parameterCountArc();

    // =============================================

    // read SST observations
    // ---------------------
    std::vector<std::vector<SatelliteTrackingArc>> sst(rhsCount);
    for(UInt j=0; j<rhsCount; j++)
    {
      sst.at(j).resize(rhs.at(j)->sstFile.size());
      for(UInt k=0; k<rhs.at(j)->sstFile.size(); k++)
      {
        sst.at(j).at(k) = rhs.at(j)->sstFile.at(k)->readArc(arcNo);
        ::Arc::checkSynchronized({sst.at(j).at(k), orbit1});
      }
    }

    // read POD observations
    // ---------------------
    std::vector<OrbitArc> pod1(rhsCount);
    std::vector<OrbitArc> pod2(rhsCount);
    for(UInt j=0; j<rhsCount; j++)
    {
      pod1.at(j) = rhs.at(j)->orbit1File->readArc(arcNo);
      pod2.at(j) = rhs.at(j)->orbit2File->readArc(arcNo);
      ::Arc::checkSynchronized({pod1.at(j), pod1.at(0)});
      ::Arc::checkSynchronized({pod2.at(j), pod2.at(0)});
    }

    // =============================================

    // count observations
    // ------------------
    const UInt countSst  = (sst.at(0).size() && (sstType!=99)) ? sst.at(0).at(0).size() : 0;
    const UInt countPod1 = 3*pod1.at(0).size();
    const UInt countPod2 = 3*pod2.at(0).size();

    UInt obsCount = 0;
    const UInt idxSst   = obsCount; obsCount += countSst;
    const UInt idxPod1  = obsCount; obsCount += countPod1;
    const UInt idxPod2  = obsCount; obsCount += countPod2;

    if(obsCount <= countBParameter)
      return Arc();

    // =============================================

    // calculate earthrotation
    // -----------------------
    std::vector<Time> times = sst.at(0).at(0).times();
    std::vector<Rotary3d> rotEarth(epochCount);
    for(UInt k=0; k<epochCount; k++)
      rotEarth.at(k) = earthRotation->rotaryMatrix(times.at(k));

    // =============================================

    // reference acceleration
    // ----------------------
    // satellite 1
    Matrix g1(3*epochCount, rhsCount);
    for(UInt j=0; j<rhsCount; j++)
    {
      AccelerometerArc accelerometer = rhs.at(j)->accelerometer1File->readArc(arcNo);
      for(UInt k=0; k<epochCount; k++)
      {
        Vector3d gv = rhs.at(j)->forces->acceleration(satellite1, orbit1.at(k).time, orbit1.at(k).position, orbit1.at(k).velocity,
                                                     starCamera1.at(k).rotary, rotEarth.at(k), earthRotation, ephemerides);
        // accelerometer [-> TRF]
        if(accelerometer.size())
          gv += rotEarth.at(k).rotate(starCamera1.at(k).rotary.rotate(accelerometer.at(k).acceleration));
        // sort into vector
        g1(3*k+0, j) = gv.x();
        g1(3*k+1, j) = gv.y();
        g1(3*k+2, j) = gv.z();
      }
    }

    // satellite 2
    Matrix g2(3*epochCount, rhsCount);
    for(UInt j=0; j<rhsCount; j++)
    {
      AccelerometerArc accelerometer = rhs.at(j)->accelerometer2File->readArc(arcNo);
      for(UInt k=0; k<epochCount; k++)
      {
        Vector3d gv = rhs.at(j)->forces->acceleration(satellite2, orbit2.at(k).time, orbit2.at(k).position, orbit2.at(k).velocity,
                                                     starCamera2.at(k).rotary, rotEarth.at(k), earthRotation, ephemerides);
        // accelerometer [-> TRF]
        if(accelerometer.size())
          gv += rotEarth.at(k).rotate(starCamera2.at(k).rotary.rotate(accelerometer.at(k).acceleration));
        // sort into vector
        g2(3*k+0, j) = gv.x();
        g2(3*k+1, j) = gv.y();
        g2(3*k+2, j) = gv.z();
      }
    }

    // =============================================

    // satellite parameters
    // --------------------
    Matrix Sat1   (3*epochCount, parameterAcceleration1->parameterCount());
    Matrix Sat2   (3*epochCount, parameterAcceleration2->parameterCount());
    Matrix Sat1Arc(3*epochCount, parameterAcceleration1->parameterCountArc());
    Matrix Sat2Arc(3*epochCount, parameterAcceleration2->parameterCountArc());
    for(UInt k=0; k<epochCount; k++)
    {
      parameterAcceleration1->compute(satellite1, times.at(k), orbit1.at(k).position, orbit1.at(k).velocity,
                                   starCamera1.at(k).rotary, rotEarth.at(k), ephemerides, Sat1.row(3*k,3), Sat1Arc.row(3*k,3));
      parameterAcceleration2->compute(satellite2, times.at(k), orbit2.at(k).position, orbit2.at(k).velocity,
                                   starCamera2.at(k).rotary, rotEarth.at(k), ephemerides, Sat2.row(3*k,3), Sat2Arc.row(3*k,3));
    }

    // =============================================

    IntegralEquation::Arc arc1 = integralEquation.integrateArc(orbit1, rotEarth, gradientfield, g1, Sat1Arc, computeVelocity || computeAcceleration, computeAcceleration);
    IntegralEquation::Arc arc2 = integralEquation.integrateArc(orbit2, rotEarth, gradientfield, g2, Sat2Arc, computeVelocity || computeAcceleration, computeAcceleration);

    // =============================================

    Matrix lSst (countSst, rhsCount);
    Matrix lSst0(countSst, rhsCount);
    Matrix VSst (countSst, 6*epochCount);
    Matrix VSstBoundary(countSst, 12);
    Matrix SstPara(countSst, parameterSst->parameterCount());
    Matrix SstParaArc(countSst, parameterSst->parameterCountArc());
    std::vector<Rotary3d> rotSat1(countSst), rotSat2(countSst);

    if(countSst)
    {
      // inter satellite vectors
      std::vector< std::vector<Vector3d> > pos12(rhsCount);
      std::vector< std::vector<Vector3d> > vel12(rhsCount);
      std::vector< std::vector<Vector3d> > acc12(rhsCount);
      std::vector< std::vector<Vector3d> > e12  (rhsCount);

      for(UInt j=0; j<rhsCount; j++)
      {
        pos12.at(j).resize(epochCount);
        vel12.at(j).resize(epochCount);
        acc12.at(j).resize(epochCount);
        e12.at(j).resize(epochCount);

        for(UInt i=0; i<epochCount; i++)
        {
          pos12.at(j).at(i) = Vector3d(arc2.vPos(3*i+0,j)-arc1.vPos(3*i+0,j),
                                       arc2.vPos(3*i+1,j)-arc1.vPos(3*i+1,j),
                                       arc2.vPos(3*i+2,j)-arc1.vPos(3*i+2,j));
          if(computeVelocity || computeAcceleration)
          vel12.at(j).at(i) = Vector3d(arc2.vVel(3*i+0,j)-arc1.vVel(3*i+0,j),
                                       arc2.vVel(3*i+1,j)-arc1.vVel(3*i+1,j),
                                       arc2.vVel(3*i+2,j)-arc1.vVel(3*i+2,j));
          if(computeAcceleration)
          acc12.at(j).at(i) = Vector3d(arc2.vAcc(3*i+0,j)-arc1.vAcc(3*i+0,j),
                                       arc2.vAcc(3*i+1,j)-arc1.vAcc(3*i+1,j),
                                       arc2.vAcc(3*i+2,j)-arc1.vAcc(3*i+2,j));

          e12.at(j).at(i) = normalize(pos12.at(j).at(i));
        }
      }

      // =============================================

      if(computeRange)
      {
        for(UInt i=0; i<epochCount; i++)
        {
          Matrix e = e12.at(0).at(i).vector().trans();
          // dRange/dPos12
          matMult(-1., e, arc1.VPos.row(3*i,3),         VSst.slice(i,0,1,3*epochCount));
          matMult( 1., e, arc2.VPos.row(3*i,3),         VSst.slice(i,3*epochCount,1,3*epochCount));
          matMult(-1., e, arc1.VPosBoundary.row(3*i,3), VSstBoundary.slice(i,0,1,6));
          matMult( 1., e, arc2.VPosBoundary.row(3*i,3), VSstBoundary.slice(i,6,1,6));
        }
      }

      if(computeVelocity)
      {
        for(UInt i=0; i<epochCount; i++)
        {
          const Double rho  = pos12.at(0).at(i).r();
          const Double drho = inner(e12.at(0).at(i),vel12.at(0).at(i));

          // dRate/dPos12
          Matrix e = ((1/rho)*vel12.at(0).at(i) - (drho/rho) * e12.at(0).at(i)).vector().trans();
          matMult(-1., e, arc1.VPos.row(3*i,3),         VSst.slice(i,0,1,3*epochCount));
          matMult( 1., e, arc2.VPos.row(3*i,3),         VSst.slice(i,3*epochCount,1,3*epochCount));
          matMult(-1., e, arc1.VPosBoundary.row(3*i,3), VSstBoundary.slice(i,0,1,6));
          matMult( 1., e, arc2.VPosBoundary.row(3*i,3), VSstBoundary.slice(i,6,1,6));

          // dRate/dVel12 = e12
          e = e12.at(0).at(i).vector().trans();
          matMult(-1., e, arc1.VVel.row(3*i,3),         VSst.slice(i,0,1,3*epochCount));
          matMult( 1., e, arc2.VVel.row(3*i,3),         VSst.slice(i,3*epochCount,1,3*epochCount));
          matMult(-1., e, arc1.VVelBoundary.row(3*i,3), VSstBoundary.slice(i,0,1,6));
          matMult( 1., e, arc2.VVelBoundary.row(3*i,3), VSstBoundary.slice(i,6,1,6));
        }
      }

      if(computeAcceleration)
      {
        for(UInt i=0; i<epochCount; i++)
        {
          const Double rho  = pos12.at(0).at(i).r();
          const Double drho = inner(e12.at(0).at(i),vel12.at(0).at(i));

          // dAccl/dPos12
          Matrix e = ((1/rho) * acc12.at(0).at(i)
                   - (2*drho/(rho*rho)) * vel12.at(0).at(i)
                   + (3*drho*drho/(rho*rho) - vel12.at(0).at(i).quadsum()/(rho*rho)-inner(e12.at(0).at(i),acc12.at(0).at(i))/rho) * e12.at(0).at(i)).vector().trans();
          matMult(-1., e, arc1.VPos.row(3*i,3),         VSst.slice(i,0,1,3*epochCount));
          matMult( 1., e, arc2.VPos.row(3*i,3),         VSst.slice(i,3*epochCount,1,3*epochCount));
          matMult(-1., e, arc1.VPosBoundary.row(3*i,3), VSstBoundary.slice(i,0,1,6));
          matMult( 1., e, arc2.VPosBoundary.row(3*i,3), VSstBoundary.slice(i,6,1,6));

          // dAccl/dVel12
          e = ((2/rho) * vel12.at(0).at(i) - (2*drho/rho) * e12.at(0).at(i)).vector().trans();
          matMult(-1., e, arc1.VVel.row(3*i,3),         VSst.slice(i,0,1,3*epochCount));
          matMult( 1., e, arc2.VVel.row(3*i,3),         VSst.slice(i,3*epochCount,1,3*epochCount));
          matMult(-1., e, arc1.VVelBoundary.row(3*i,3), VSstBoundary.slice(i,0,1,6));
          matMult( 1., e, arc2.VVelBoundary.row(3*i,3), VSstBoundary.slice(i,6,1,6));

          // dAccl/dAcc12 = e12
          e = e12.at(0).at(i).vector().trans();
          matMult(-1., e, arc1.VAcc.row(3*i,3),         VSst.slice(i,0,1,3*epochCount));
          matMult( 1., e, arc2.VAcc.row(3*i,3),         VSst.slice(i,3*epochCount,1,3*epochCount));
          matMult(-1., e, arc1.VAccBoundary.row(3*i,3), VSstBoundary.slice(i,0,1,6));
          matMult( 1., e, arc2.VAccBoundary.row(3*i,3), VSstBoundary.slice(i,6,1,6));
        }
      }

      // =============================================

      // reference observations
      // ----------------------
      for(UInt j=0; j<rhsCount; j++)
        for(UInt i=0; i<epochCount; i++)
        {
          const Double rho  = pos12.at(j).at(i).r();
          const Double drho = inner(e12.at(j).at(i), vel12.at(j).at(i));

          if(computeRange)        lSst0(i,j) = rho;
          if(computeVelocity)     lSst0(i,j) = drho;
          if(computeAcceleration) lSst0(i,j) = inner(e12.at(j).at(i),acc12.at(j).at(i)) + (1/rho) * (vel12.at(j).at(i).quadsum() - drho*drho);
        }

      // reduced observations
      // ----------------------
      for(UInt j=0; j<rhsCount; j++)
        for(UInt k=0; k<sst.at(j).size(); k++)
          for(UInt i=0; i<epochCount; i++)
          {
            if(computeRange)        lSst(i,j) += sst.at(j).at(k).at(i).range;
            if(computeVelocity)     lSst(i,j) += sst.at(j).at(k).at(i).rangeRate;
            if(computeAcceleration) lSst(i,j) += sst.at(j).at(k).at(i).rangeAcceleration;
          }
      lSst -= lSst0;

      // =============================================

      // satellite tracking parameter
      // ----------------------------
      for(UInt i=0; i<countSst; i++)
        rotSat1.at(i) = starCamera1.at(i).rotary;
      for(UInt i=0; i<countSst; i++)
        rotSat2.at(i) = starCamera2.at(i).rotary;

      parameterSst->compute(sstType, times, lSst0.column(0), arc1.vPos, arc2.vPos, arc1.vVel, arc2.vVel, rotSat1, rotSat2, SstPara, SstParaArc);
    } // if(countSst)


    // =============================================

    Matrix lPos1, VPos1, VPosBoundary1;
    Matrix lPos2, VPos2, VPosBoundary2;
    integralEquation.interpolateArc(pod1, orbit1, arc1, lPos1, VPos1, VPosBoundary1);
    integralEquation.interpolateArc(pod2, orbit2, arc2, lPos2, VPos2, VPosBoundary2);
    arc1 = IntegralEquation::Arc();
    arc2 = IntegralEquation::Arc();

    // =============================================

    // observation vector
    // ------------------
    Matrix l(obsCount, rhsCount);
    if(countSst)  copy(lSst,  l.row(idxSst,  countSst));
    if(countPod1) copy(lPos1, l.row(idxPod1, countPod1));
    if(countPod2) copy(lPos2, l.row(idxPod2, countPod2));

    // =============================================

    // integration of accelerations to observations
    // --------------------------------------------
    Matrix V(obsCount, 6*epochCount);
    if(countSst)  copy(VSst,  V.slice(idxSst,  0,            countSst,  6*epochCount));
    if(countPod1) copy(VPos1, V.slice(idxPod1, 0,            countPod1, 3*epochCount));
    if(countPod2) copy(VPos2, V.slice(idxPod2, 3*epochCount, countPod2, 3*epochCount));
    VSst = VPos1 = VPos2 = Matrix();

    // =============================================

    // arc related parameters (satellite state + acceleration calibration)
    // -------------------------------------------------------------------
    Matrix B(obsCount, countBParameter); // 4 boundary pos. (x,y,z) + AccCal + CalKBandArc

    // satellite parameters
    if(Sat1Arc.size()) matMult(1., V.column(0,           3*epochCount), Sat1Arc, B.column(idxSat1Arc, Sat1Arc.columns()));
    if(Sat2Arc.size()) matMult(1., V.column(3*epochCount,3*epochCount), Sat2Arc, B.column(idxSat2Arc, Sat2Arc.columns()));
    Sat1Arc = Sat2Arc = Matrix();

    // satellite tracking parameters
    if(SstParaArc.size()) copy(SstParaArc, B.slice(idxSst, idxSstParaArc, countSst, SstParaArc.columns()));
    SstParaArc = Matrix();

    // boundaries (satellite state vector)
    if(!keepSatelliteStates)
    {
      if(countSst)  copy(VSstBoundary.column(0,6), B.slice(idxSst,  idxBound1, countSst,  6));
      if(countSst)  copy(VSstBoundary.column(6,6), B.slice(idxSst,  idxBound2, countSst,  6));
      if(countPod1) copy(VPosBoundary1,            B.slice(idxPod1, idxBound1, countPod1, 6));
      if(countPod2) copy(VPosBoundary2,            B.slice(idxPod2, idxBound2, countPod2, 6));
      VSstBoundary = VPosBoundary1 = VPosBoundary2 = Matrix();
    }

    // =============================================

    // Design matrix A (calibration part)
    // ----------------------------------
    Matrix A(obsCount, countAParameter);

    // satellite parameters
    if(Sat1.size()) matMult(1., V.column(0,           3*epochCount), Sat1, A.column(idxSat1, Sat1.columns()));
    if(Sat2.size()) matMult(1., V.column(3*epochCount,3*epochCount), Sat2, A.column(idxSat2, Sat2.columns()));
    Sat1 = Sat2 = Matrix();

    // satellite tracking parameters
    if(SstPara.size()) copy(SstPara, A.slice(idxSst, idxSstPara, countSst, SstPara.columns()));
    SstPara = Matrix();

    // boundaries (satellite state vector)
    if(keepSatelliteStates)
    {
      if(countSst)  copy(VSstBoundary.column(0,6), A.slice(idxSst,  idxBound1+6*arcNo, countSst,  6));
      if(countSst)  copy(VSstBoundary.column(6,6), A.slice(idxSst,  idxBound2+6*arcNo, countSst,  6));
      if(countPod1) copy(VPosBoundary1,            A.slice(idxPod1, idxBound1+6*arcNo, countPod1, 6));
      if(countPod2) copy(VPosBoundary2,            A.slice(idxPod2, idxBound2+6*arcNo, countPod2, 6));
      VSstBoundary = VPosBoundary1 = VPosBoundary2 = Matrix();
    }

    // =============================================

    // decorrelation
    // -------------
    if(covSst  && countSst)  covSst->decorrelate (arcNo, times,      {l.row(idxSst,  countSst),  V.row(idxSst,  countSst),  B.row(idxSst,  countSst),  A.slice(idxSst,  gravityCount, countSst,  countAParameter-gravityCount)});
    if(covPod1 && countPod1) covPod1->decorrelate(arcNo, pod1.at(0), {l.row(idxPod1, countPod1), V.row(idxPod1, countPod1), B.row(idxPod1, countPod1), A.slice(idxPod1, gravityCount, countPod1, countAParameter-gravityCount)});
    if(covPod2 && countPod2) covPod2->decorrelate(arcNo, pod2.at(0), {l.row(idxPod2, countPod2), V.row(idxPod2, countPod2), B.row(idxPod2, countPod2), A.slice(idxPod2, gravityCount, countPod2, countAParameter-gravityCount)});

    // =============================================

    // Design matrix A (gravity)
    // -------------------------
    if(gravityCount)
    {
      Matrix G(3*epochCount, gravityCount);
      for(UInt k=0; k<epochCount; k++)
        parameterGravity->gravity(orbit1.at(k).time, rotEarth.at(k).rotate(orbit1.at(k).position), G.row(3*k,3));
      matMult(1., V.column(0,3*epochCount), G, A.column(idxGravity, gravityCount));

      for(UInt k=0; k<epochCount; k++)
        parameterGravity->gravity(orbit2.at(k).time, rotEarth.at(k).rotate(orbit2.at(k).position), G.row(3*k,3));
      matMult(1., V.column(3*epochCount,3*epochCount), G, A.column(idxGravity, gravityCount));
    }

    // =============================================

    Arc observationArc;
    observationArc.l     = l;
    observationArc.A     = A;
    observationArc.B     = B;
    observationArc.times = {times, pod1.front().times(), pod2.front().times()};
    observationArc.pod1  = pod1.front();
    observationArc.pod2  = pod2.front();
    return observationArc;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
