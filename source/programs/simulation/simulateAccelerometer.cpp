/***********************************************/
/**
* @file simulateAccelerometer.cpp
*
* @brief simulate accelerometer data.
*
* @author Torsten Mayer-Guerr
* @date 2005-01-21
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program simulate \file{accelerometer data}{instrument}. The orientation of the accelerometer
is given by \configFile{inputfileStarCamera}{instrument} otherwise the celestial reference frame (CRF) is used.
For computation of non-conservative forces a \configFile{satelliteModel}{satelliteModel} is needed.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "files/fileSatelliteModel.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/ephemerides/ephemerides.h"
#include "classes/forces/forces.h"

/***** CLASS ***********************************/

/** @brief Simulate accelerometer data.
* @ingroup programsGroup */
class SimulateAccelerometer
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(SimulateAccelerometer, PARALLEL, "simulate accelerometer data", Simulation, Instrument)

/***********************************************/

void SimulateAccelerometer::run(Config &config)
{
  try
  {
    FileName             fileNameAccelerometer;
    FileName             fileNameSatellite;
    FileName             orbitName, starCameraName;
    EarthRotationPtr     earthRotation;
    EphemeridesPtr       ephemerides;
    ForcesPtr            forces;

    renameDeprecatedConfig(config, "satelliteModel", "inputfileSatelliteModel", date2time(2020, 8, 19));

    readConfig(config, "outputfileAccelerometer", fileNameAccelerometer, Config::MUSTSET,  "",    "");
    readConfig(config, "inputfileSatelliteModel", fileNameSatellite,     Config::OPTIONAL, "{groopsDataDir}/satelliteModel/", "satellite macro model");
    readConfig(config, "inputfileOrbit",          orbitName,             Config::MUSTSET,  "",    "");
    readConfig(config, "inputfileStarCamera",     starCameraName,        Config::OPTIONAL, "",    "");
    readConfig(config, "earthRotation",           earthRotation,         Config::MUSTSET,  "",    "");
    readConfig(config, "ephemerides",             ephemerides,           Config::OPTIONAL, "jpl", "");
    readConfig(config, "forces",                  forces,                Config::MUSTSET,  "",    "");
    if(isCreateSchema(config)) return;

    InstrumentFile orbitFile(orbitName);
    InstrumentFile starCameraFile(starCameraName);
    InstrumentFile::checkArcCount({orbitFile, starCameraFile});

    SatelliteModelPtr satellite;
    if(!fileNameSatellite.empty())
      readFileSatelliteModel(fileNameSatellite, satellite);

    logStatus<<"computing accelerations"<<Log::endl;
    std::vector<Arc> arcList(orbitFile.arcCount());
    Parallel::forEach(arcList, [&](UInt arcNo)
    {
      OrbitArc      orbit      = orbitFile.readArc(arcNo);
      StarCameraArc starCamera = starCameraFile.readArc(arcNo);
      Arc::checkSynchronized({orbit, starCamera});

      AccelerometerArc accelerometer;
      for(UInt k=0; k<orbit.size(); k++)
      {
        Rotary3d rotSat;
        if(starCamera.size())
          rotSat = starCamera.at(k).rotary;
        const Time     time     = orbit.at(k).time;
        const Rotary3d rotEarth = earthRotation->rotaryMatrix(time);
        const Vector3d acc      = forces->acceleration(satellite, time, orbit.at(k).position, orbit.at(k).velocity, rotSat, rotEarth, earthRotation, ephemerides);

        AccelerometerEpoch epoch;
        epoch.time         = time;
        epoch.acceleration = rotSat.inverseRotate(rotEarth.inverseRotate(acc));
        accelerometer.push_back(epoch);
      }
      return accelerometer;
    });

    if(Parallel::isMaster())
    {
      logStatus<<"write accelerometer data to file <"<<fileNameAccelerometer<<">"<<Log::endl;
      InstrumentFile::write(fileNameAccelerometer, arcList);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
