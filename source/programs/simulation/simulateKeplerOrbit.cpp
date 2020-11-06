/***********************************************/
/**
* @file simulateKeplerOrbit.cpp
*
* @brief Compute Keplerian orbit.
*
* @author Torsten Mayer-Guerr
* @date 2009-10-31
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program simulates a Keplerian \file{orbit}{instrument} at a given \config{timeSeries}
starting from the given \config{integrationConstants}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/kepler.h"
#include "base/equinoctial.h"
#include "files/fileInstrument.h"
#include "classes/timeSeries/timeSeries.h"

/***** CLASS ***********************************/

/** @brief Compute Keplerian orbit.
* @ingroup programsGroup */
class SimulateKeplerOrbit
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(SimulateKeplerOrbit, SINGLEPROCESS, "compute Keplerian orbit", Simulation, Orbit, Instrument)

/***********************************************/

void SimulateKeplerOrbit::run(Config &config)
{
  try
  {
    FileName      outName;
    Double        GM;
    TimeSeriesPtr timeSeries;
    Vector3d      position, velocity;
    std::string   choice;
    Time          time0;

    readConfig(config, "outputfileOrbit", outName,       Config::MUSTSET,  "", "");
    readConfig(config, "timeSeries",      timeSeries,    Config::MUSTSET,  "", "");
    readConfig(config, "GM",              GM,            Config::DEFAULT,  STRING_DEFAULT_GM, "Geocentric gravitational constant");
    readConfigChoice(config, "integrationConstants", choice, Config::MUSTSET, "", "");
    if(readConfigChoiceElement(config, "kepler", choice, ""))
    {
      Double a, e;
      Angle  i, Omega, omega, M;
      readConfig(config, "majorAxis",         a,     Config::MUSTSET, "", "[m]");
      readConfig(config, "eccentricity",      e,     Config::MUSTSET, "", "[-]");
      readConfig(config, "inclination",       i,     Config::MUSTSET, "", "[degree]");
      readConfig(config, "ascendingNode",     Omega, Config::MUSTSET, "", "[degree]");
      readConfig(config, "argumentOfPerigee", omega, Config::MUSTSET, "", "[degree]");
      readConfig(config, "meanAnomaly",       M,     Config::MUSTSET, "", "[degree]");
      readConfig(config, "time",              time0, Config::MUSTSET, "", "integration constants are valid at this epoch");
      if(!isCreateSchema(config))
      {
        Kepler kepler(time0, Omega, i, omega, a, e, M, GM);
        kepler.orbit(time0, position, velocity);
      }
    }
    if(readConfigChoiceElement(config, "positionAndVelocity", choice, ""))
    {
      readConfig(config, "position0x", position.x(), Config::MUSTSET, "", "[m] in CRF");
      readConfig(config, "position0y", position.y(), Config::MUSTSET, "", "[m] in CRF");
      readConfig(config, "position0z", position.z(), Config::MUSTSET, "", "[m] in CRF");
      readConfig(config, "velocity0x", velocity.x(), Config::MUSTSET, "", "[m/s]");
      readConfig(config, "velocity0y", velocity.y(), Config::MUSTSET, "", "[m/s]");
      readConfig(config, "velocity0z", velocity.z(), Config::MUSTSET, "", "[m/s]");
      readConfig(config, "time",       time0,        Config::MUSTSET, "", "integration constants are valid at this epoch");
    }
    endChoice(config);
    if(isCreateSchema(config)) return;

    std::vector<Time> times = timeSeries->times();
    Equinoctial kepler(time0, position, velocity, GM);

    // Info
    // ----
    Kepler k(kepler);
    logInfo<<"Kepler elements"<<Log::endl;
    logInfo<<"  majorAxis a    : "<<k.a/1000<<" km"<<Log::endl;
    logInfo<<"  eccentricity e : "<<k.e<<Log::endl;
    logInfo<<"  inclination i  : "<<k.i*RAD2DEG<<" Degree"<<Log::endl;
    logInfo<<"  ascending node : "<<k.Omega*RAD2DEG<<" Degree"<<Log::endl;
    logInfo<<"  perigee        : "<<k.omega*RAD2DEG<<" Degree"<<Log::endl;

    // Computation
    // -----------
    logStatus<<"computing"<<Log::endl;
    OrbitArc orbit;
    logTimerStart;
    for(UInt i=0; i<times.size(); i++)
    {
      logTimerLoop(i,times.size());
      OrbitEpoch epoch;
      epoch.time = times.at(i);
      kepler.orbit(times.at(i), epoch.position, epoch.velocity, epoch.acceleration);
      orbit.push_back(epoch);
    }
    logTimerLoopEnd(times.size());

    // write results
    // -------------
    logStatus<<"write orbit data to file <"<<outName<<">"<<Log::endl;
    InstrumentFile::write(outName, orbit);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
