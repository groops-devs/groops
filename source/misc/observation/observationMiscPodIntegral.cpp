/***********************************************/
/**
* @file observationMiscPodIntegral.cpp
*
* @brief Precise Orbit data (Short Arc Integral).
*
* @author Torsten Mayer-Guerr
* @date 2015-06-02
*
*/
/***********************************************/

#include "base/import.h"
#include "parallel/parallel.h"
#include "files/fileInstrument.h"
#include "files/fileSatelliteModel.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/gravityfield/gravityfield.h"
#include "classes/parametrizationGravity/parametrizationGravity.h"
#include "classes/parametrizationAcceleration/parametrizationAcceleration.h"
#include "misc/observation/observationMisc.h"
#include "observationMiscPodIntegral.h"

/***********************************************/

ObservationMiscPodIntegral::ObservationMiscPodIntegral(Config &config)
{
  try
  {
    FileName fileNameSatellite;
    FileName orbitName, starCameraName;

    renameDeprecatedConfig(config, "satelliteModel", "inputfileSatelliteModel",     date2time(2020, 8, 19));
    renameDeprecatedConfig(config, "representation", "parametrizationGravity",      date2time(2020, 6,  3));
    renameDeprecatedConfig(config, "parameter",      "parametrizationAcceleration", date2time(2020, 6,  3));

    readConfig(config, "inputfileSatelliteModel",     fileNameSatellite,     Config::OPTIONAL, "{groopsDataDir}/satelliteModel/", "satellite macro model");
    readConfig(config, "rightHandSide",               rhs,                   Config::MUSTSET,  "",    "input for the reduced observation vector");
    readConfig(config, "inputfileOrbit",              orbitName,             Config::MUSTSET,  "",    "used to evaluate the observation equations, not used as observations");
    readConfig(config, "inputfileStarCamera",         starCameraName,        Config::MUSTSET,  "",    "");
    readConfig(config, "earthRotation",               earthRotation,         Config::MUSTSET,  "",    "");
    readConfig(config, "ephemerides",                 ephemerides,           Config::OPTIONAL, "jpl", "");
    readConfig(config, "gradientfield",               gradientfield,         Config::DEFAULT,  "",    "low order field to estimate the change of the gravity by position adjustement");
    readConfig(config, "parametrizationGravity",      parameterGravity,      Config::DEFAULT,  "",    "gravity field parametrization");
    readConfig(config, "parametrizationAcceleration", parameterAcceleration, Config::DEFAULT,  "",    "orbit force parameters");
    readConfig(config, "keepSatelliteStates",         keepSatelliteStates,   Config::DEFAULT,  "0",   "set boundary values of each arc global");
    readConfig(config, "integrationDegree",           integrationDegree,     Config::DEFAULT,  "7",   "integration of forces by polynomial approximation of degree n");
    readConfig(config, "interpolationDegree",         interpolationDegree,   Config::DEFAULT,  "7",   "orbit interpolation by polynomial approximation of degree n");
    readConfig(config, "accelerateComputation",       accelerateComputation, Config::DEFAULT,  "0",   "acceleration of computation by transforming the observations");
    if(isCreateSchema(config)) return;

    if(integrationDegree%2 == 0)
      throw(Exception("polnomial degree for integration must be odd."));

    if(!fileNameSatellite.empty())
      readFileSatelliteModel(fileNameSatellite, satellite);

    // test instrument files
    // ---------------------
    orbitFile.open(orbitName);
    starCameraFile.open(starCameraName);
    InstrumentFile::checkArcCount({orbitFile, starCameraFile});
    for(UInt j=0; j<rhs.size(); j++)
      InstrumentFile::checkArcCount({orbitFile, *(rhs.at(j)->accelerometerFile)});
    for(UInt j=1; j<rhs.size(); j++)
      InstrumentFile::checkArcCount({*(rhs.at(j)->orbitFile), *(rhs.at(0)->orbitFile)});

    // init
    // ----
    countArc = orbitFile.arcCount();
    integralEquation.init(integrationDegree, interpolationDegree);

    // design matrix A
    // ---------------
    countAParameter = 0;
    idxGravity = countAParameter; countAParameter += parameterGravity->parameterCount();
    idxSat     = countAParameter; countAParameter += parameterAcceleration->parameterCount();
    if(keepSatelliteStates)
    {
      idxBound = countAParameter;
      countAParameter += 6*countArc; // 2 boundary pos. (x,y,z).
    }

    if(accelerateComputation && keepSatelliteStates)
    {
      accelerateComputation = FALSE;
      if(Parallel::isMaster())
        logWarning<<"acceleration of computation is not possible, if keepSatelliteStates is set"<<Log::endl;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool ObservationMiscPodIntegral::setInterval(const Time &timeStart, const Time &timeEnd)
{
  try
  {
    Bool change = FALSE;
    change = parameterGravity->setInterval(timeStart, timeEnd)      || change;
    change = parameterAcceleration->setInterval(timeStart, timeEnd) || change;

    // count parameters
    // ----------------
    countAParameter = 0;
    idxGravity = countAParameter; countAParameter += parameterGravity->parameterCount();
    idxSat     = countAParameter; countAParameter += parameterAcceleration->parameterCount();
    if(keepSatelliteStates)
    {
      idxBound = countAParameter;
      countAParameter += 6*countArc; // 2 boundary pos. (x,y,z).
    }

    return change;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ObservationMiscPodIntegral::parameterName(std::vector<ParameterName> &name) const
{
  try
  {
    parameterGravity->parameterName(name);

    const std::string satelliteName = satellite ? satellite->satelliteName : "satellite";
    parameterAcceleration->parameterName(name);
    for(UInt i = name.size(); i --> name.size()-parameterAcceleration->parameterCount(); )
      name.at(i).object = satelliteName;

    if(keepSatelliteStates)
    {
      for(UInt arcNo=0; arcNo<arcCount(); arcNo++)
      {
        const std::string str = "arc"+arcNo%"%i"s+".position.";
        name.push_back(ParameterName(satelliteName, str+"start.x"));
        name.push_back(ParameterName(satelliteName, str+"start.y"));
        name.push_back(ParameterName(satelliteName, str+"start.z"));
        name.push_back(ParameterName(satelliteName, str+"end.x"));
        name.push_back(ParameterName(satelliteName, str+"end.y"));
        name.push_back(ParameterName(satelliteName, str+"end.z"));
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

ObservationMiscPod::Arc ObservationMiscPodIntegral::computeArc(UInt arcNo, CovariancePodPtr covPod)
{
  try
  {
    OrbitArc      orbit      = orbitFile.readArc(arcNo);
    StarCameraArc starCamera = starCameraFile.readArc(arcNo);
    const UInt    rhsCount   = rhs.size();
    const UInt    epochCount = orbit.size();

    parameterAcceleration->setIntervalArc(orbit.at(0).time, orbit.back().time+medianSampling(orbit.times()));

    // count arc related parameters (in matrix B)
    // ------------------------------------------
    UInt countBParameter = 0;
    if(!keepSatelliteStates)
    {
      idxBound = countBParameter;
      countBParameter += 6;  // 2 Randpos. je (x,y,z).
    }
    const UInt idxSatArc = countBParameter;
    countBParameter += parameterAcceleration->parameterCountArc();

    // read POD observations
    // ---------------------
    std::vector<OrbitArc> pod(rhsCount);
    for(UInt j=0; j<rhsCount; j++)
      pod.at(j) = rhs.at(j)->orbitFile->readArc(arcNo);

    for(UInt j=1; j<rhsCount; j++)
      for(UInt i=0; i<pod.at(j).size(); i++)
        if(pod.at(j).at(i).time != pod.at(j-1).at(i).time)
          throw(Exception("POD-Orbits on different rightHandSides differ in sampling"));

    const UInt podCount = pod.at(0).size();
    if(3*podCount <= countBParameter)
      return Arc();

    // calculate earthrotation
    // -----------------------
    std::vector<Rotary3d> rotEarth(epochCount);
    for(UInt k=0; k<epochCount; k++)
      rotEarth.at(k) = earthRotation->rotaryMatrix(orbit.at(k).time);

    // reference acceleration
    // ----------------------
    Matrix g(3*epochCount, rhsCount);
    AccelerometerArc accelerometer;
    for(UInt j=0; j<rhsCount; j++)
    {
      AccelerometerArc accl = rhs.at(j)->accelerometerFile->readArc(arcNo);
      if((accelerometer.size() == 0) && (accl.size() != 0))
        accelerometer = accl;

      for(UInt k=0; k<epochCount; k++)
      {
        Vector3d gv = rhs.at(j)->forces->acceleration(satellite, orbit.at(k).time, orbit.at(k).position, orbit.at(k).velocity,
                                                      starCamera.at(k).rotary, rotEarth.at(k), earthRotation, ephemerides);
        // accelerometer [-> TRF]
        if(accl.size())
          gv += rotEarth.at(k).rotate(starCamera.at(k).rotary.rotate(accl.at(k).acceleration));
        // sort into vector
        g(3*k+0, j) = gv.x();
        g(3*k+1, j) = gv.y();
        g(3*k+2, j) = gv.z();
      }
    }

    // =============================================

    // satellite parameters
    // --------------------
    Matrix Sat   (3*epochCount, parameterAcceleration->parameterCount());
    Matrix SatArc(3*epochCount, parameterAcceleration->parameterCountArc());
    for(UInt k=0; k<epochCount; k++)
    {
      parameterAcceleration->compute(satellite, orbit.at(k).time, orbit.at(k).position, orbit.at(k).velocity,
                                  starCamera.at(k).rotary, rotEarth.at(k), ephemerides, Sat.row(3*k,3), SatArc.row(3*k,3));
    }

    // =============================================

    // integral kernel
    // ---------------
    Matrix l;
    Matrix VPos, VPosBoundary;
    {
      IntegralEquation::Arc arc = integralEquation.integrateArc(orbit, rotEarth, gradientfield, g, SatArc);
      integralEquation.interpolateArc(pod, orbit, arc, l, VPos, VPosBoundary);
    }

    // =============================================

    // decorrelation
    // -------------
    if(covPod)
      covPod->decorrelate(arcNo, pod.at(0), {l, VPos, VPosBoundary});

    // =============================================

    // arc related parameters (satellite state + acceleration calibration)
    // -------------------------------------------------------------------
    Matrix B(3*podCount, countBParameter);
    // satellite state vector
    if(!keepSatelliteStates)
      copy(VPosBoundary, B.column(idxBound, VPosBoundary.columns()));
    // Accelerometer calibration
    if(SatArc.columns() != 0)
      matMult(1., VPos, SatArc, B.column(idxSatArc, SatArc.columns()));

    // =============================================

    if(accelerateComputation)
    {
      eliminationParameter(B, VPos, l);
      B = Matrix();
      if(VPos.rows() > VPos.columns())
      {
        Vector tau = QR_decomposition(VPos);
        QTransMult(VPos, tau, l);
        VPos = VPos.row(0,VPos.columns());
        VPos.setType(Matrix::TRIANGULAR, Matrix::UPPER);
      }
    }

    // =============================================

    // Design matrix A (gravity)
    // -------------------------
    Matrix A(VPos.rows(), countAParameter);
    // gravity field
    if(parameterGravity->parameterCount())
    {
      Matrix G(3*epochCount, parameterGravity->parameterCount());
      for(UInt k=0; k<epochCount; k++)
        parameterGravity->gravity(orbit.at(k).time, rotEarth.at(k).rotate(orbit.at(k).position), G.row(3*k, 3));
      matMult(1., VPos, G, A.column(idxGravity, G.columns()));
    }
    // Accelerometer calibration
    if(Sat.columns() != 0)
      matMult(1., VPos, Sat, A.column(idxSat, Sat.columns()));
    // satellite state vector
    if(keepSatelliteStates)
      copy(VPosBoundary, A.column(idxBound+6*arcNo, VPosBoundary.columns()));

    // =============================================

    Arc observationArc;
    observationArc.l     = l;
    observationArc.A     = A;
    observationArc.B     = B;
    observationArc.times = pod.at(0).times();
    observationArc.pod   = pod.at(0);
    return observationArc;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
