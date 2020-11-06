/***********************************************/
/**
* @file observationPodAcceleration.cpp
*
* @brief Precise Orbit Data (POD) observations (Acceleration approach).
*
* @author Torsten Mayer-Guerr
* @date 2006-02-02
* update 2011-06-09 Norbert Zehentner
*
*/
/***********************************************/

#include "base/import.h"
#include "parallel/parallel.h"
#include "files/fileInstrument.h"
#include "files/fileSatelliteModel.h"
#include "misc/observation/observationMisc.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/parametrizationGravity/parametrizationGravity.h"
#include "classes/parametrizationAcceleration/parametrizationAcceleration.h"
#include "classes/observation/observation.h"
#include "classes/observation/observationPodAcceleration.h"

/***********************************************/

ObservationPodAcceleration::ObservationPodAcceleration(Config &config)
{
  try
  {
    FileName fileNameSatellite;
    FileName orbitName, starCameraName;
    UInt     epochCount, interpolationDegree;

    renameDeprecatedConfig(config, "satelliteModel", "inputfileSatelliteModel",     date2time(2020, 8, 19));
    renameDeprecatedConfig(config, "representation", "parametrizationGravity",      date2time(2020, 6,  3));
    renameDeprecatedConfig(config, "parameter",      "parametrizationAcceleration", date2time(2020, 6,  3));

    readConfig(config, "inputfileSatelliteModel",     fileNameSatellite,            Config::OPTIONAL, "{groopsDataDir}/satelliteModel/", "satellite macro model");
    readConfig(config, "rightHandSide",               rhs,                          Config::MUSTSET,  "",    "input for the reduced observation vector");
    readConfig(config, "inputfileOrbit",              orbitName,                    Config::MUSTSET,  "",    "used to evaluate the observation equations, not used as observations");
    readConfig(config, "inputfileStarCamera",         starCameraName,               Config::MUSTSET,  "",    "");
    readConfig(config, "earthRotation",               earthRotation,                Config::MUSTSET,  "",    "");
    readConfig(config, "ephemerides",                 ephemerides,                  Config::OPTIONAL, "jpl", "");
    readConfig(config, "parametrizationGravity",      parametrization,              Config::MUSTSET,  "",    "gravity field parametrization");
    readConfig(config, "parametrizationAcceleration", parametrizationAcceleration,  Config::DEFAULT,  "",    "orbit force parameters");
    readConfig(config, "interpolationDegree",         interpolationDegree,          Config::DEFAULT,  "8",   "orbit differentation  by polynomial approximation of degree n");
    readConfig(config, "numberOfEpochs",              epochCount,                   Config::DEFAULT,  "9",   "number of used Epochs for polynom computation");
    readConfig(config, "covariancePod",               covPod,                       Config::OPTIONAL, "",    "covariance matrix of kinematic orbits");
    if(isCreateSchema(config)) return;

    // test if epochCount and polynomial degree are correct
    // ---------------------
    if(interpolationDegree%2 != 0)
      throw(Exception("polnomial degree for interpolation must be even."));
    if(epochCount <= interpolationDegree)
      throw(Exception("Epochs must be greater than interpolationdegree, at least by one."));
    if(epochCount%2 != 1)
      throw(Exception("Epochs must be odd (symmetry around central point)."));

    if(!fileNameSatellite.empty())
      readFileSatelliteModel(fileNameSatellite, satellite);

    // test instrument files
    // ---------------------
    orbitFile.open(orbitName);
    starCameraFile.open(starCameraName);
    InstrumentFile::checkArcCount({orbitFile, starCameraFile});
    for(UInt j=0; j<rhs.size(); j++)
      InstrumentFile::checkArcCount({orbitFile, *rhs.at(j)->orbitFile, *rhs.at(j)->accelerometerFile});

    // orbit differentation coefficients (position -> acceleration) by means of QR-decomposition
    // with or without overdetermination
    // -----------------------
    Matrix P(epochCount, interpolationDegree+1);
    for(UInt i=0; i<epochCount; i++)
      for(UInt n=0; n<=interpolationDegree; n++)
        P(i,n) = ((n==0) ? 1.0 : std::pow(i-(epochCount-1.)/2., n));

    const Vector tau = QR_decomposition(P);
    coeff = Vector(epochCount);
    coeff(2) = 2.0;
    triangularSolve(1., P.row(0,interpolationDegree+1).trans(), coeff);
    QMult(P, tau, coeff);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

UInt ObservationPodAcceleration::parameterCount() const
{
  return parametrization->parameterCount() + parametrizationAcceleration->parameterCount();
}

/***********************************************/

void ObservationPodAcceleration::parameterName(std::vector<ParameterName> &name) const
{
  if(parametrization)
    parametrization->parameterName(name);

  if(parametrizationAcceleration)
  {
    parametrizationAcceleration->parameterName(name);
    const std::string satelliteName = satellite ? satellite->satelliteName : "satellite";
    for(UInt i=name.size(); i-->name.size()-parametrizationAcceleration->parameterCount();)
      name.at(i).object = satelliteName;
  }
}

/***********************************************/

void ObservationPodAcceleration::observation(UInt arcNo, Matrix &l, Matrix &A, Matrix &B)
{
  try
  {
    OrbitArc      orbit      = orbitFile.readArc(arcNo);
    StarCameraArc starCamera = starCameraFile.readArc(arcNo);
    UInt          rhsCount   = rhs.size();
    UInt          posCount   = orbit.size();
    UInt          obsCount   = posCount-coeff.rows()+1;
    UInt          half       = (coeff.rows()-1)/2;
    Double        dt         = medianSampling(orbit.times()).seconds();
    Arc::checkSynchronized({orbit, starCamera});

    parametrizationAcceleration->setIntervalArc(orbit.at(half).time, orbit.at(obsCount+half-1).time+seconds2time(dt));

    // calculate earthrotation
    // -----------------------
    std::vector<Rotary3d> rotEarth(obsCount);
    for(UInt k=0; k<obsCount; k++)
      rotEarth.at(k) = earthRotation->rotaryMatrix(orbit.at(k+half).time);

    // reference acceleration
    // ----------------------
    std::vector<std::vector<Vector3d>> reference(rhsCount);
    for(UInt j=0; j<rhsCount; j++)
    {
      reference.at(j).resize(obsCount);
      AccelerometerArc accelerometer = rhs.at(j)->accelerometerFile->readArc(arcNo);
      Arc::checkSynchronized({orbit, accelerometer});

      for(UInt k=0; k<obsCount; k++)
      {
        Vector3d g = rhs.at(j)->forces->acceleration(satellite, orbit.at(k+half).time, orbit.at(k+half).position, orbit.at(k+half).velocity,
                                                     starCamera.at(k+half).rotary, rotEarth.at(k), earthRotation, ephemerides);
        if(accelerometer.size())
          g += rotEarth.at(k).rotate(starCamera.at(k+half).rotary.rotate(accelerometer.at(k+half).acceleration)); // accelerometer
        reference.at(j).at(k) = g; // reference observations
      }
    }

    // reduced observations
    // --------------------
    l = Matrix(3*obsCount, rhsCount);
    for(UInt j=0; j<rhsCount; j++)
    {
      OrbitArc orbitPod = rhs.at(j)->orbitFile->readArc(arcNo);
      for(UInt k=0; k<obsCount; k++)
        {
          Vector3d a;
          for(UInt i=0; i<coeff.rows(); i++)
            a += (coeff(i)/dt/dt) * orbitPod.at(k+i).position;
          a = rotEarth.at(k).rotate(a);

          l(3*k+0,j) = a.x() - reference.at(j).at(k).x();
          l(3*k+1,j) = a.y() - reference.at(j).at(k).y();
          l(3*k+2,j) = a.z() - reference.at(j).at(k).z();
        }
    }

    // Design matrix A and B
    // ---------------------
    A = Matrix(3*obsCount, parametrization->parameterCount() + parametrizationAcceleration->parameterCount());
    B = Matrix(3*obsCount, parametrizationAcceleration->parameterCountArc());
    // gravity
    for(UInt k=0; k<obsCount; k++)
      parametrization->gravity(orbit.at(k+half).time, rotEarth.at(k).rotate(orbit.at(k+half).position), A.slice(3*k, 0, 3, parametrization->parameterCount()));

    // orbit parameters
    for(UInt k=0; k<obsCount; k++)
    {
      parametrizationAcceleration->compute(satellite, orbit.at(k+half).time, orbit.at(k+half).position, orbit.at(k+half).velocity,
                                           starCamera.at(k+half).rotary, rotEarth.at(k), ephemerides,
                                           A.slice(3*k, parametrization->parameterCount(), 3, parametrizationAcceleration->parameterCount()), B.row(3*k,3));
    }

    // decorrelation
    // -------------
    if(covPod)
    {
      // linearized variance propagation
      Matrix D(3*obsCount, 3*posCount);
      for(UInt k=0; k<obsCount; k++)
        for(UInt i=0; i<coeff.rows(); i++)
          axpy((coeff(i)/dt/dt), rotEarth.at(k).matrix(), D.slice(3*k,3*(k+i),3,3));
      Matrix C = covPod->covariance(arcNo, orbit);
      Matrix DCD = D * C * D.trans();
      DCD.setType(Matrix::SYMMETRIC);
      cholesky(DCD);

      // apply Cholesky matrix
      triangularSolve(1., DCD.trans(), A);
      triangularSolve(1., DCD.trans(), l);
      if(B.size())
        triangularSolve(1., DCD.trans(), B);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
