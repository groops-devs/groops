/***********************************************/
/**
* @file observationPodEnergy.cpp
*
* @brief Precise Orbit Data (POD) observations (Energy integral).
*
* @author Torsten Mayer-Guerr
* @date 2006-02-02
*
*/
/***********************************************/

#include "base/import.h"
#include "parallel/parallel.h"
#include "files/fileInstrument.h"
#include "files/fileSatelliteModel.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/parametrizationGravity/parametrizationGravity.h"
#include "classes/observation/observation.h"
#include "misc/observation/observationMisc.h"
#include "misc/observation/covariancePod.h"
#include "classes/observation/observationPodEnergy.h"

/***********************************************/

ObservationPodEnergy::ObservationPodEnergy(Config &config)
{
  try
  {
    FileName fileNameSatellite;
    FileName orbitName, starCameraName;

    renameDeprecatedConfig(config, "satelliteModel", "inputfileSatelliteModel", date2time(2020, 8, 19));
    renameDeprecatedConfig(config, "representation", "parametrizationGravity",  date2time(2020, 6,  3));

    readConfig(config, "inputfileSatelliteModel", fileNameSatellite,   Config::OPTIONAL, "{groopsDataDir}/satelliteModel/", "satellite macro model");
    readConfig(config, "rightHandSide",           rhs,                 Config::MUSTSET,   "",   "input for the reduced observation vector");
    readConfig(config, "inputfileOrbit",          orbitName,           Config::MUSTSET,   "",   "used to evaluate the observation equations, not used as observations");
    readConfig(config, "inputfileStarCamera",     starCameraName,      Config::MUSTSET,   "",   "");
    readConfig(config, "earthRotation",           earthRotation,       Config::MUSTSET,   "",   "");
    readConfig(config, "ephemerides",             ephemerides,         Config::OPTIONAL, "jpl", "");
    readConfig(config, "parametrizationGravity",  parametrization,     Config::MUSTSET,   "",   "gravity field parametrization (potential)");
    readConfig(config, "parametrizationBias",     bias,                Config::MUSTSET,   "",   "unknown total energy per arc");
    readConfig(config, "interpolationDegree",     interpolationDegree, Config::DEFAULT,   "8",  "orbit differentation  by polynomial approximation of degree n");
    readConfig(config, "integrationDegree",       integrationDegree,   Config::DEFAULT,   "7",  "integration of forces by polynomial approximation of degree n");
    readConfig(config, "covariancePod",           covPod,              Config::OPTIONAL,  "",   "covariance matrix of kinematic orbits");
    if(isCreateSchema(config)) return;

    if(interpolationDegree%2 != 0)
      throw(Exception("polnomial degree for interpolation must be even."));

    if(!fileNameSatellite.empty())
      readFileSatelliteModel(fileNameSatellite, satellite);

    // test instrument files
    // ---------------------
    orbitFile.open(orbitName);
    starCameraFile.open(starCameraName);
    InstrumentFile::checkArcCount({orbitFile, starCameraFile});
    for(UInt j=0; j<rhs.size(); j++)
      InstrumentFile::checkArcCount({orbitFile, *rhs.at(j)->orbitFile, *rhs.at(j)->accelerometerFile});

    // orbit differentation coefficients (position -> velocity)
    Matrix P(interpolationDegree+1, interpolationDegree+1);
    for(UInt i=0; i<=interpolationDegree; i++)
      for(UInt n=0; n<=interpolationDegree; n++)
        P(n,i) = ((n==0) ? 1.0 : std::pow((static_cast<Double>(i)-interpolationDegree/2), n));
    coeff = Vector(interpolationDegree+1);
    coeff(1) = 1.0;
    solveInPlace(P, coeff);

    // polynomial integration matrix
    integrationMatrix = Matrix(integrationDegree+1, integrationDegree+1);
    for(UInt i=0; i<integrationMatrix.rows(); i++)
    {
      integrationMatrix(0,i) = 1.0;
      for(UInt n=1; n<integrationMatrix.columns(); n++)
        integrationMatrix(n,i) = (i-integrationDegree/2.) * integrationMatrix(n-1,i);
    }
    inverse(integrationMatrix);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ObservationPodEnergy::observation(UInt arcNo, Matrix &l, Matrix &A, Matrix &B)
{
  try
  {
    OrbitArc      orbit      = orbitFile.readArc(arcNo);
    StarCameraArc starCamera = starCameraFile.readArc(arcNo);
    UInt          rhsCount   = rhs.size();
    UInt          posCount   = orbit.size();
    UInt          obsCount   = posCount-interpolationDegree;
    UInt          half       = interpolationDegree/2;
    Double        dt         = (orbit.at(half+1).time - orbit.at(half).time).seconds();
    Arc::checkSynchronized({orbit, starCamera});

    // calculate earthrotation
    // -----------------------
    std::vector<Rotary3d>  rotEarth(obsCount);
    std::vector<Vector3d>  rotAxis(obsCount);
    std::vector<Vector3d>  rotAxisDot(obsCount);
    for(UInt k=0; k<obsCount; k++)
    {
      const Time time = orbit.at(k+half).time;
      rotEarth.at(k)   = earthRotation->rotaryMatrix(time);
      rotAxis.at(k)    = earthRotation->rotaryAxis(time);
      rotAxisDot.at(k) = earthRotation->rotaryAxisDerivate(time);
    }

    // position + velocity
    // -------------------
    std::vector<std::vector<Vector3d>> position(rhsCount);
    std::vector<std::vector<Vector3d>> velocity(rhsCount);
    std::vector<std::vector<Vector3d>> velocityEarth(rhsCount);
    for(UInt j=0; j<rhsCount; j++)
    {
      OrbitArc orbit = rhs.at(j)->orbitFile->readArc(arcNo);
      position.at(j).resize(obsCount);
      velocity.at(j).resize(obsCount);
      velocityEarth.at(j).resize(obsCount);
      for(UInt k=0; k<obsCount; k++)
      {
        position.at(j).at(k) = orbit.at(k+half).position;
        for(UInt i=0; i<coeff.rows(); i++)
          velocity.at(j).at(k) += (coeff(i)/dt) * orbit.at(k+i).position;
        // Earth fixed velocity in spaced fixed coordinates
        velocityEarth.at(j).at(k) = velocity.at(j).at(k) - crossProduct(rotAxis.at(k), position.at(j).at(k));
      }
    }

    // reference acceleration
    // ----------------------
    Matrix integrand(obsCount, rhsCount);
    for(UInt j=0; j<rhsCount; j++)
    {
      AccelerometerArc accl = rhs.at(j)->accelerometerFile->readArc(arcNo);
      for(UInt k=0; k<obsCount; k++)
      {
        Vector3d g = rhs.at(j)->forces->acceleration(satellite, orbit.at(k+half).time, orbit.at(k+half).position, orbit.at(k+half).velocity,
                                                     starCamera.at(k+half).rotary, rotEarth.at(k), earthRotation, ephemerides);
        g = rotEarth.at(k).inverseRotate(g);                               // rotation into CRF
        if(accl.size()!=0)
          g += starCamera.at(k+half).rotary.rotate(accl.at(k+half).acceleration); // accelerometer

        integrand(k,j) = inner(g, velocityEarth.at(j).at(k))
                       - inner(crossProduct(rotAxisDot.at(k), position.at(j).at(k)), velocity.at(j).at(k));
      }
    }

    // reduced observations
    // --------------------
    l = Matrix(obsCount, rhsCount);
    for(UInt j=0; j<rhsCount; j++)
    {
      const Double l0 = 0.5 * velocity.at(j).at(0).quadsum() - inner(velocity.at(j).at(0), crossProduct(rotAxis.at(0), position.at(j).at(0)));
      for(UInt k=1; k<obsCount; k++)
        l(k,j) = 0.5 * velocity.at(j).at(k).quadsum() - inner(velocity.at(j).at(k), crossProduct(rotAxis.at(k), position.at(j).at(k))) - l0;
    }

    // reduce integral
    // ---------------
    Matrix integral(integrand.rows(), integrand.columns());
    Vector integralFactors;
    for(UInt i=1; i<integrand.rows(); i++)
    {
      const UInt idx = std::min(std::max(i, integrationDegree/2)-integrationDegree/2, integrand.rows()-integrationDegree-1);
      if((i<=integrationDegree/2) || (i>=integrand.rows()-(integrationDegree+1)/2))
      {
        integralFactors = Vector(integrationDegree+1);
        const Double tau1 = i-idx-integrationDegree/2.-1;
        const Double tau2 = i-idx-integrationDegree/2.;
        for(UInt n=0; n<integralFactors.rows(); n++)
          axpy(dt/(n+1)*(std::pow(tau2, n+1)-std::pow(tau1, n+1)), integrationMatrix.column(n), integralFactors);
      }
      integral.row(i) += integral.row(i-1);
      matMult(-1., integralFactors.trans(), integrand.row(idx, integralFactors.rows()), integral.row(i));
    }
    l -= integral;

    // arc related parameters (energy constant)
    // ----------------------------------------
    B = Matrix();
    const std::vector<Time> times = orbit.times();
    bias->setInterval(times.front(), times.back()+medianSampling(times), TRUE);
    if(bias->parameterCount())
    {
      B = Matrix(obsCount, bias->parameterCount());
      for(UInt i=0; i<obsCount; i++)
        copy(bias->factors(orbit.at(i+half).time).trans(), B.row(i));
    }

    // Design matrix A (gravitational potential)
    // -----------------------------------------
    A = Matrix(obsCount, parametrization->parameterCount());
    for(UInt k=0; k<obsCount; k++)
      parametrization->potential(orbit.at(k+half).time, rotEarth.at(k).rotate(orbit.at(k+half).position), A.row(k));

    // decorrelation
    // -------------
    if(covPod)
    {
      // linearized variance propagation
      Matrix D(obsCount, 3*posCount);
      for(UInt k=0; k<obsCount; k++)
        for(UInt i=0; i<coeff.rows(); i++)
        {
          D(k,3*(i+k)+0) = (coeff(i)/dt) * velocity.at(0).at(k).x();
          D(k,3*(i+k)+1) = (coeff(i)/dt) * velocity.at(0).at(k).y();
          D(k,3*(i+k)+2) = (coeff(i)/dt) * velocity.at(0).at(k).z();
        }

      Matrix C = covPod->covariance(arcNo, orbit);
      Matrix DCD = D * C * D.trans();
      DCD.setType(Matrix::SYMMETRIC);
      cholesky(DCD);

      // apply Cholesky matrix
      triangularSolve(1., DCD.trans(), A);
      triangularSolve(1., DCD.trans(), l);
      if(B.size()!=0)
        triangularSolve(1., DCD.trans(), B);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
