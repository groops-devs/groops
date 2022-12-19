/***********************************************/
/**
* @file instrumentAccelerometer2ThermosphericDensity.cpp
*
* @brief Estimate neutral density from accelerometer data.
*
* @author Sandro Krauss
* @date 2020-06-22
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program estimates neutral mass densities along the satellite trajectory based on \file{accelerometer data}{instrument}.
In order to determine the neutral mass density the accelerometer input should only reflect the accelerations due to drag
(e.g. \configClass{miscAccelerations:atmosphericDrag}{miscAccelerationsType:atmosphericDrag}).
Thus, influences from solar and Earth radiation pressure must be reduced beforehand.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileMatrix.h"
#include "files/fileSatelliteModel.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/parametrizationAcceleration/parametrizationAcceleration.h"
#include "classes/thermosphere/thermosphere.h"
#include "classes/miscAccelerations/miscAccelerationsAtmosphericDrag.h"
#include "classes/ephemerides/ephemerides.h"

/***** CLASS ***********************************/

/** @brief Estimate neutral density from accelerometer data.
* @ingroup programGroup */
class InstrumentAccelerometer2ThermosphericDensity
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(InstrumentAccelerometer2ThermosphericDensity, PARALLEL, "Estimate neutral density from accelerometer data", Instrument)

/***********************************************/

void InstrumentAccelerometer2ThermosphericDensity::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName         fileNameOutDensity;
    FileName         fileNameSatellite;
    FileName         fileNameOrbit, fileNameStarCamera, fileNameAccelerometer;
    EarthRotationPtr earthRotation;
    ThermospherePtr  thermosphere;
    Bool             useTemperature, useWind;
    EphemeridesPtr   ephemerides;

    readConfig(config, "outputfileDensity",      fileNameOutDensity,     Config::MUSTSET,  "",    "MISCVALUE (kg/m^3)");
    readConfig(config, "satelliteModel",         fileNameSatellite,      Config::OPTIONAL, "{groopsDataDir}/satelliteModel/", "satellite macro model");
    readConfig(config, "inputfileOrbit",         fileNameOrbit,          Config::MUSTSET,  "",    "");
    readConfig(config, "inputfileStarCamera",    fileNameStarCamera,     Config::MUSTSET,  "",    "");
    readConfig(config, "inputfileAccelerometer", fileNameAccelerometer,  Config::MUSTSET,  "",    "add non-gravitational forces in satellite reference frame");
    readConfig(config, "thermosphere",           thermosphere,           Config::MUSTSET,  "",    "used to compute temperature and wind");
    readConfig(config, "considerTemperature",    useTemperature,         Config::DEFAULT,  "1",   "compute drag and lift, otherwise simple drag coefficient is used");
    readConfig(config, "considerWind",           useWind,                Config::DEFAULT,  "1",   "");
    readConfig(config, "earthRotation",          earthRotation,          Config::MUSTSET,  "",    "");
    readConfig(config, "ephemerides",            ephemerides,            Config::OPTIONAL, "jpl", "");
    if(isCreateSchema(config)) return;

    // open and test instrument files
    // ------------------------------
    InstrumentFile orbitFile, starCameraFile, accelerometerFile;
    orbitFile.open(fileNameOrbit);
    starCameraFile.open(fileNameStarCamera);
    accelerometerFile.open(fileNameAccelerometer);
    InstrumentFile::checkArcCount({orbitFile, starCameraFile, accelerometerFile});

    SatelliteModelPtr satellite;
    if(!fileNameSatellite.empty())
      readFileSatelliteModel(fileNameSatellite, satellite);

    logStatus<<"computing density"<<Log::endl;
    std::vector<Arc> arcList(orbitFile.arcCount());
    Parallel::forEach(arcList, [&](UInt arcNo)
    {
      OrbitArc      orbit      = orbitFile.readArc(arcNo);
      StarCameraArc starCamera = starCameraFile.readArc(arcNo);
      AccelerometerArc accelerometer = accelerometerFile.readArc(arcNo);
      Arc::checkSynchronized({orbit, starCamera, accelerometer});

      MiscValueArc output;
      for(UInt k=0; k<orbit.size(); k++)
      {
        Rotary3d rotSat;
        if(starCamera.size())
          rotSat = starCamera.at(k).rotary;
        const Time     time     = orbit.at(k).time;
        const Rotary3d rotEarth = earthRotation->rotaryMatrix(time);
        const Vector3d omega    = earthRotation->rotaryAxis(orbit.at(k).time);

        if(satellite)
        {
          Vector3d positionSun;
          if(ephemerides)
            positionSun = ephemerides->position(time, Ephemerides::SUN);
          satellite->changeState(time, orbit.at(k).position, orbit.at(k).velocity, positionSun, rotSat, rotEarth);
        }

        Double   modelDensity, temperature;
        Vector3d wind;
        thermosphere->state(time, rotEarth.rotate(orbit.at(k).position), modelDensity, temperature, wind);
        if(!useTemperature)
         temperature = 0;
        if(!useWind)
         wind = Vector3d();

        // direction and speed of thermosphere relative to satellite in SRF
        Vector3d direction = rotSat.inverseRotate(rotEarth.inverseRotate(wind) + crossProduct(omega, orbit.at(k).position) - orbit.at(k).velocity);
        const Double v = direction.normalize();
        const Vector3d acc = (1./satellite->mass) * MiscAccelerationsAtmosphericDrag::force(satellite, direction, v, 1., temperature);

        MiscValueEpoch epoch;
        epoch.time  = time;
        epoch.value = accelerometer.at(k).acceleration.x() / acc.x(); // density
        output.push_back(epoch);
      }
      return output;
    }, comm);

    if(Parallel::isMaster(comm))
    {
      logStatus<<"write density data to file <"<<fileNameOutDensity<<">"<<Log::endl;
      InstrumentFile::write(fileNameOutDensity, arcList);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
