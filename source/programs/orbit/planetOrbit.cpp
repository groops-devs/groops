/***********************************************/
/**
* @file planetOrbit.cpp
*
* @brief Orbits of sun, moon, and planets.
*
* @author Torsten Mayer-Guerr
* @date 2011-04-12
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Creates an \file{orbit file}{instrument} of sun, moon, or planets.
The orbit is given in the celestial reference frame (CRF)
or alternatively in the terrestrial refernce frame (TRF)
if \configClass{earthRotation}{earthRotationType} is provided.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "classes/timeSeries/timeSeries.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/ephemerides/ephemerides.h"

/***** CLASS ***********************************/

/** @brief Orbits of sun, moon, and planets.
* @ingroup programsGroup */
class PlanetOrbit
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(PlanetOrbit, SINGLEPROCESS, "orbits of sun, moon and, planets", Orbit, Instrument)

/***********************************************/

void PlanetOrbit::run(Config &config)
{
  try
  {
    FileName            fileNameOrbit;
    TimeSeriesPtr       timeSeries;
    EarthRotationPtr    earthRotation;
    EphemeridesPtr      ephemerides;
    Ephemerides::Planet planet;
    std::string         choice;

    readConfig(config, "outputfileOrbit", fileNameOrbit, Config::MUSTSET,  "", "");
    readConfig(config, "planet",          planet,        Config::MUSTSET,  "", "");
    readConfig(config, "timeSeries",      timeSeries,    Config::MUSTSET,  "", "");
    readConfig(config, "ephemerides",     ephemerides,   Config::MUSTSET,  "", "");
    readConfig(config, "earthRotation",   earthRotation, Config::OPTIONAL, "", "transform orbits into TRF");
    if(isCreateSchema(config)) return;

    // =============================================

    logStatus<<"computing"<<Log::endl;
    std::vector<Time> times  = timeSeries->times();
    OrbitArc orbit;
    logTimerStart;
    for(UInt i=0; i<times.size(); i++)
    {
      logTimerLoop(i,times.size());
      OrbitEpoch epoch;
      epoch.time = times.at(i);
      ephemerides->ephemeris(times.at(i), planet, epoch.position, epoch.velocity);

      if(earthRotation)
      {
        const Rotary3d rotEarth = earthRotation->rotaryMatrix(times.at(i));
        const Vector3d omega    = earthRotation->rotaryAxis(times.at(i));
        epoch.velocity = rotEarth.rotate(epoch.velocity - crossProduct(omega, epoch.position));
        epoch.position = rotEarth.rotate(epoch.position);
      }

      orbit.push_back(epoch);
    }
    logTimerLoopEnd(times.size());

    logStatus<<"write orbit data to file <"<<fileNameOrbit<<">"<<Log::endl;
    InstrumentFile::write(fileNameOrbit, orbit);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
