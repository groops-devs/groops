/***********************************************/
/**
* @file simulateAccelerometerCoMOffset.cpp
*
* @brief Create accelerations due to CoM offset.
*
* @author Beate Klinger
* @date 2016-05-12
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program generates an \file{accelerometer file}{instrument} containing perturbing accelerations
due to a given center of mass (CoM) offset. This includes centrifugal effects,
Euler forces and the effect of gravity gradients.
)";

/***********************************************/

#include "programs/program.h"
#include "base/polynomial.h"
#include "files/fileInstrument.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/gravityfield/gravityfield.h"

/***** CLASS ***********************************/

/** @brief Create accelerations due to CoM offset.
* @ingroup programsGroup */
class SimulateAccelerometerCoMOffset
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(SimulateAccelerometerCoMOffset, PARALLEL, "Create accelerations due to CoM offset.", Simulation, Instrument)

/***********************************************/

void SimulateAccelerometerCoMOffset::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName         fileNameOut;
    FileName         fileNameOrbit, fileNameStarCamera;
    Bool             applyAngularRate, applyAngularAcc;
    GravityfieldPtr  gradientfield;
    EarthRotationPtr earthRotation;
    UInt             derivationDegree;
    Vector           offset(3);

    readConfig(config, "outputfileAccelerometer",   fileNameOut,        Config::MUSTSET,  "", "effect of offset");
    readConfig(config, "inputfileOrbit",            fileNameOrbit,      Config::MUSTSET,  "", "");
    readConfig(config, "inputfileStarCamera",       fileNameStarCamera, Config::MUSTSET,  "", "");
    readConfig(config, "applyAngularRate",          applyAngularRate,   Config::DEFAULT,  "1", "compute effect of centrifugal forces");
    readConfig(config, "applyAngularAccelerations", applyAngularAcc,    Config::DEFAULT,  "1", "compute effect of Euler forces");
    readConfig(config, "gradientfield",             gradientfield,      Config::DEFAULT,  "", "low order field to estimate the change of the gravity by position adjustement");
    readConfig(config, "earthRotation",             earthRotation,      Config::MUSTSET,  "", "");
    readConfig(config, "interpolationDegree",       derivationDegree,   Config::DEFAULT,  "2", "derivation of quaternions by polynomial interpolation of degree n");
    readConfig(config, "CoMOffsetX",                offset(0),          Config::MUSTSET,  "50e-6", "offset [m]");
    readConfig(config, "CoMOffsetY",                offset(1),          Config::MUSTSET,  "50e-6", "offset [m]");
    readConfig(config, "CoMOffsetZ",                offset(2),          Config::MUSTSET,  "50e-6", "offset [m]");
    if(isCreateSchema(config)) return;

    Polynomial polynomial(derivationDegree);

    logStatus<<"compute accelerations (CoM)"<<Log::endl;
    InstrumentFile orbitFile(fileNameOrbit);
    InstrumentFile starCameraFile(fileNameStarCamera);
    InstrumentFile::checkArcCount({orbitFile, starCameraFile});
    std::vector<Arc> arcList(orbitFile.arcCount());

    Parallel::forEach(arcList, [&](UInt arcNo)
    {
      OrbitArc         orbit      = orbitFile.readArc(arcNo);
      StarCameraArc    starCamera = starCameraFile.readArc(arcNo);
      Arc::checkSynchronized({orbit, starCamera});

      // quaternions
      Double dt  = medianSampling(orbit.times()).seconds();
      Matrix q   = starCamera.matrix().column(1,4);
      Matrix dq  = polynomial.derivative(dt, q);
      Matrix ddq = polynomial.derivative2nd(dt, q);

      AccelerometerArc accelerometerArc;
      for(UInt i=0; i<orbit.size(); i++)
      {
        // Gravity gradient in satellite system (SRF)
        Rotary3d rotEarth = earthRotation->rotaryMatrix(orbit.at(i).time);
        Tensor3d T = gradientfield->gravityGradient(orbit.at(i).time, rotEarth.rotate(orbit.at(i).position));
        T = starCamera.at(i).rotary.inverseRotate(rotEarth.inverseRotate(T));

        // angular rates/accelerations
        Vector w(3), wd(3);
        Matrix Q(3, 4);
        Q(0,0) = -q(i,1);  Q(0,1) = +q(i,0);  Q(0,2) = +q(i,3);  Q(0,3) = -q(i,2);
        Q(1,1) = -q(i,2);  Q(1,1) = -q(i,3);  Q(1,2) = +q(i,0);  Q(1,3) = +q(i,1);
        Q(2,2) = -q(i,3);  Q(2,1) = +q(i,2);  Q(2,2) = -q(i,1);  Q(2,3) = +q(i,0);
        if(applyAngularRate)
          w = 2*Q *  dq.row(i).trans();
        if(applyAngularAcc)
          wd = 2*Q * ddq.row(i).trans(); // angular accelerations

        // Center of mass offset
        Matrix M = -T.matrix();
        M(0,0) += -w(1)*w(1)-w(2)*w(2); M(0,1) +=  w(0)*w(1)-wd(2);     M(0,2) +=  w(0)*w(2)+wd(1);
        M(1,0) +=  w(0)*w(1)+wd(2);     M(1,1) += -w(0)*w(0)-w(2)*w(2); M(1,2) +=  w(1)*w(2)-wd(0);
        M(2,0) +=  w(0)*w(2)-wd(1);     M(2,1) +=  w(1)*w(2)+wd(0);     M(2,2) += -w(0)*w(0)-w(1)*w(1);
        Vector a = M * offset;

        AccelerometerEpoch epoch;
        epoch.time = orbit.at(i).time;
        epoch.acceleration.x() = a(0);
        epoch.acceleration.y() = a(1);
        epoch.acceleration.z() = a(2);
        accelerometerArc.push_back(epoch);
      }
      return accelerometerArc;
    }, comm);

    if(Parallel::isMaster(comm))
    {
      logStatus<<"write data to file <"<<fileNameOut<<">"<<Log::endl;
      InstrumentFile::write(fileNameOut, arcList);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
