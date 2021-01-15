/***********************************************/
/**
* @file simulateGradiometer.cpp
*
* @brief Simulate error free gradiometer data.
*
* @author Torsten Mayer-Guerr
* @date 2008-08-15
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program simulates error free \file{gradiometer data}{instrument} along a satellite's orbit.
The orientation of the full tensor gradiometer is given by \configFile{inputfileStarCamera}{instrument}
otherwise the celestial reference frame (CRF) is used.
The gravity gradients are given by \configClass{gravityfield}{gravityfieldType} and
\configClass{tides}{tidesType}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/ephemerides/ephemerides.h"
#include "classes/tides/tides.h"
#include "classes/gravityfield/gravityfield.h"

/***********************************************/

/** @brief Simulate error free gradiometer data.
* @ingroup programsGroup */
class SimulateGradiometer
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(SimulateGradiometer, PARALLEL, "simulate error free gradiometer data", Simulation, Instrument)

/***********************************************/

void SimulateGradiometer::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName         fileNameGradiometer;
    FileName         orbitName, starCameraName;
    EarthRotationPtr earthRotation;
    EphemeridesPtr   ephemerides;
    GravityfieldPtr  gravityfield;
    TidesPtr         tides;

    readConfig(config, "outputfileGradiometer", fileNameGradiometer, Config::MUSTSET,  "", "");
    readConfig(config, "inputfileOrbit",        orbitName,           Config::MUSTSET,  "", "");
    readConfig(config, "inputfileStarCamera",   starCameraName,      Config::OPTIONAL, "", "");
    readConfig(config, "earthRotation",         earthRotation,       Config::MUSTSET,  "", "");
    readConfig(config, "ephemerides",           ephemerides,         Config::OPTIONAL, "jpl", "");
    readConfig(config, "gravityfield",          gravityfield,        Config::DEFAULT,  "", "");
    readConfig(config, "tides",                 tides,               Config::DEFAULT,  "", "");
    if(isCreateSchema(config)) return;

    // open and test instrument files
    // ------------------------------
    InstrumentFile orbitFile(orbitName);
    InstrumentFile starCameraFile(starCameraName);
    InstrumentFile::checkArcCount({orbitFile, starCameraFile});

    logStatus<<"computing gravity field"<<Log::endl;
    std::vector<Arc> arcList(orbitFile.arcCount());
    Parallel::forEach(arcList, [&](UInt arcNo)
    {
      OrbitArc      orbit      = orbitFile.readArc(arcNo);
      StarCameraArc starCamera = starCameraFile.readArc(arcNo);
      Arc::checkSynchronized({orbit, starCamera});

      GradiometerArc gradiometer;
      for(UInt k=0; k<orbit.size(); k++)
      {
        Rotary3d rotSat;
        if(starCamera.size())
          rotSat = starCamera.at(k).rotary;
        const Time     time     = orbit.at(k).time;
        const Rotary3d rotEarth = earthRotation->rotaryMatrix(time);
        const Vector3d posEarth = rotEarth.rotate(orbit.at(k).position);
        const Tensor3d tns      = gravityfield->gravityGradient(time, posEarth)
                                + tides->gradient(time, posEarth, rotEarth, earthRotation, ephemerides);

        GradiometerEpoch epoch;
        epoch.time            = time;
        epoch.gravityGradient = rotSat.inverseRotate(rotEarth.inverseRotate(tns));
        gradiometer.push_back(epoch);
      }
      return gradiometer;
    }, comm);

    if(Parallel::isMaster(comm))
    {
      logStatus<<"write gradiometer file <"<<fileNameGradiometer<<">"<<Log::endl;
      InstrumentFile::write(fileNameGradiometer, arcList);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
