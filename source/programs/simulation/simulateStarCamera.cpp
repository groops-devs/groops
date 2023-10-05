/***********************************************/
/**
* @file simulateStarCamera.cpp
*
* @brief Rotation from satellite (x: along, y: cross, z: nadir) to inertial frame.
*
* @author Torsten Mayer-Guerr
* @date 2005-01-21
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program simulates \file{star camera}{instrument} measurements at each satellite's position.
The satellite's orientation follows a local orbit frame with the x-axis in along track (along velocity),
y-axis is cross track (normal to position and velocity vector) and z-axis pointing nadir (negative position vector).
As for non circular orbit the position and velocity are not exact normal, the default is the x-axis to be exact
along velocity and the z-axis forms a right hand system (not exact nadir) or with \config{nadirPointing} the z-axis
is exact nadir and x-axis approximates along.
The resulting rotation matrices rotate from satellite frame to inertial frame.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Rotation from satellite (x: along, y: cross, z: nadir) to inertial frame.
* @ingroup programsGroup */
class SimulateStarCamera
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(SimulateStarCamera, PARALLEL, "Rotation from satellite (x: along, y: cross, z: nadir) to inertial frame.", Simulation, Instrument)

/***********************************************/

void SimulateStarCamera::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName orbitName, starCameraName;
    Bool     isNadirPointing;

    readConfig(config, "outputfileStarCamera", starCameraName,  Config::MUSTSET, "",  "rotation from satellite to inertial frame (x: along, y: cross, z: nadir)");
    readConfig(config, "inputfileOrbit",       orbitName,       Config::MUSTSET, "",  "position and velocity defines the orientation of the satellite at each epoch");
    readConfig(config, "nadirPointing",        isNadirPointing, Config::DEFAULT, "0", "false: exact along and nearly nadir, true: nearly along and exact nadir");
    if(isCreateSchema(config)) return;

    logStatus<<"read orbit and generate star camera data"<<Log::endl;
    InstrumentFile  orbitFile(orbitName);
    std::vector<Arc> arcList(orbitFile.arcCount());

    Parallel::forEach(arcList, [&](UInt arcNo)
    {
      OrbitArc orbit = orbitFile.readArc(arcNo);
      StarCameraArc arc;
      for(UInt i=0; i<orbit.size(); i++)
      {

        Vector3d velocity = orbit.at(i).velocity; // velocity vector
        if(velocity.r()==0)
        {
          if(i<orbit.size()-1)
            velocity = orbit.at(i+1).position - orbit.at(i).position;
          else
            velocity = orbit.at(i).position - orbit.at(i-1).position;
        }

        Vector3d y = normalize(crossProduct(velocity, orbit.at(i).position)); // cross
        Vector3d x = normalize((isNadirPointing) ? crossProduct(orbit.at(i).position, y) : velocity);  // along

        StarCameraEpoch epoch;
        epoch.time   = orbit.at(i).time;
        epoch.rotary = Rotary3d(x, y);
        arc.push_back(epoch);
      }
      return arc;
    }, comm);

    if(Parallel::isMaster(comm))
    {
      logStatus<<"write star camera data to file <"<<starCameraName<<">"<<Log::endl;
      InstrumentFile::write(starCameraName, arcList);
      Arc::printStatistics(arcList);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
