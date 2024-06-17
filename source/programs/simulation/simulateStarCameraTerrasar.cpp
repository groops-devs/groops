/***********************************************/
/**
* @file simulateStarCameraTerrasar.cpp
*
* @brief Simulate star camera data for Terrasar.
*
* @author Torsten Mayer-Guerr
* @date 2024-06-04
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program simulates \configFile{outputfileStarCamera}{instrument} measurements at each satellite's position for the Terrasar satellite.
The \configFile{inputfileOrbit}{instrument} must contain positions and velocities (see \program{OrbitAddVelocityAndAcceleration}).
The resulting rotation matrices rotate from satellite frame to inertial frame.

H. Fiedler, E. Boerner, J. Mittermayer and G. Krieger,
Total zero Doppler Steering-a new method for minimizing the Doppler centroid,
in IEEE Geoscience and Remote Sensing Letters, vol. 2, no. 2, pp. 141-145, April 2005, \url{https://www.doi.org/10.1109/LGRS.2005.844591}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/kepler.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief simulate star camera data for Terrasar.
* @ingroup programsGroup */
class SimulateStarCameraTerrasar
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(SimulateStarCameraTerrasar, PARALLEL, "simulate star camera data for Terrasar", Simulation, Instrument)

/***********************************************/

void SimulateStarCameraTerrasar::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName fileNameOrbit, fileNameStarCamera;

    readConfig(config, "outputfileStarCamera", fileNameStarCamera, Config::MUSTSET, "", "rotation from satellite to inertial frame (x: along, y: cross, z: nadir)");
    readConfig(config, "inputfileOrbit",       fileNameOrbit,      Config::MUSTSET, "", "position and velocity defines the orientation of the satellite at each epoch");
    if(isCreateSchema(config)) return;

    // StarCamera-Daten erzeugen
    // -----------------------
    logStatus<<"read orbit and generate star camera data"<<Log::endl;
    InstrumentFile orbitFile(fileNameOrbit);
    std::vector<Arc> arcList(orbitFile.arcCount());
    Parallel::forEach(arcList, [&](UInt arcNo)
    {
      OrbitArc orbit = orbitFile.readArc(arcNo);
      StarCameraArc arc;
      for(UInt i=0; i<orbit.size(); i++)
      {
        // Argument of latitude of satellite
        const Vector3d z = normalize(crossProduct(orbit.at(i).position, orbit.at(i).velocity));
        const Vector3d x = normalize(crossProduct(Vector3d(0,0,1), z));
        const Vector3d y = crossProduct(z, x);
        const Double   u = atan2(inner(orbit.at(i).position, y), inner(orbit.at(i).position, x));

        Kepler k(orbit.at(i).time, orbit.at(i).position, orbit.at(i).velocity);
        const Double N   = 167./11.; // 86400/std::sqrt(std::pow(k.a, 3)/k.GM);  // number of revolutions per day
        const Double yaw = std::atan((std::sin(k.i)*std::cos(u))/(N-std::cos(k.i)));

        StarCameraEpoch epoch;
        epoch.time   = orbit.at(i).time;
        epoch.rotary = Rotary3d(normalize(orbit.at(i).velocity), -z) * rotaryZ(Angle(yaw));
        arc.push_back(epoch);
      }
      return arc;
    }, comm); // forEach

    if(Parallel::isMaster(comm))
    {
      logStatus<<"write star camera data to file <"<<fileNameStarCamera<<">"<<Log::endl;
      InstrumentFile::write(fileNameStarCamera, arcList);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
