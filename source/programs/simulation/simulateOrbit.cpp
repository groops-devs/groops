/***********************************************/
/**
* @file simulateOrbit.cpp
*
* @brief Pure dynamical orbit integration.
*
* @author Matthias Ellmer
* @date 2004-09-12
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program integrates an \file{orbit}{instrument} from a given force function (dynamic orbit).
The force functions are given by \configClass{forces}{forcesType}.
For computation of non-conservative forces a \file{satelliteModel}{satelliteModel} is needed.
The integration method must be selected with \configClass{propagator}{orbitPropagatorType}.
Because the orbit data are calculated in the celestial reference frame (CRF) you need
\configClass{earthRotation}{earthRotationType} to transform the force function
from the terrestrial reference frame (TRF).
The integration start and end time, as well as the sampling, are derived from
the \config{timeSeries} option. It is possible to integrate the arc in \config{reverse},
where the initial conditions are assumed to be met at the end time of the \config{timeSeries}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/kepler.h"
#include "files/fileInstrument.h"
#include "files/fileSatelliteModel.h"
#include "classes/timeSeries/timeSeries.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/ephemerides/ephemerides.h"
#include "classes/forces/forces.h"
#include "classes/orbitPropagator/orbitPropagator.h"

/***** CLASS ***********************************/

/** @brief Pure dynamical orbit integration.
 * @ingroup programsGroup */
class SimulateOrbit
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(SimulateOrbit, SINGLEPROCESS, "pure dynamical orbit integration", Simulation, Orbit, Instrument)

/***********************************************/

void SimulateOrbit::run(Config &config)
{
  try
  {
    FileName             fileNameOrbit;
    FileName             fileNameSatellite;
    FileName             fileNameStart;
    TimeSeriesPtr        timeSeries;
    OrbitPropagatorPtr   orbitPropagator;
    OrbitEpoch           startEpoch;
    EarthRotationPtr     earthRotation;
    EphemeridesPtr       ephemerides;
    ForcesPtr            forces;
    Bool                 reverse;
    Double               margin;

    renameDeprecatedConfig(config, "satelliteModel", "inputfileSatelliteModel", date2time(2020, 8, 19));

    readConfig(config, "outputfileOrbit",         fileNameOrbit,     Config::MUSTSET,   "",   "orbit file to be written.");
    readConfig(config, "inputfileSatelliteModel", fileNameSatellite, Config::OPTIONAL, "{groopsDataDir}/satelliteModel/", "satellite macro model");
    readConfig(config, "timeSeries",              timeSeries,        Config::MUSTSET,   "",   "time points for simulated orbit epochs.");

    std::string choice;
    readConfigChoice(config, "integrationConstants", choice, Config::MUSTSET, "", "");
    if(readConfigChoiceElement(config, "kepler", choice, ""))
    {
      Double GM, a, e;
      Angle  i, Omega, omega, M;
      readConfig(config, "majorAxis",         a,     Config::MUSTSET, "", "[m]");
      readConfig(config, "eccentricity",      e,     Config::MUSTSET, "", "[-]");
      readConfig(config, "inclination",       i,     Config::MUSTSET, "", "[degree]");
      readConfig(config, "ascendingNode",     Omega, Config::MUSTSET, "", "[degree]");
      readConfig(config, "argumentOfPerigee", omega, Config::MUSTSET, "", "[degree]");
      readConfig(config, "meanAnomaly",       M,     Config::MUSTSET, "", "[degree]");
      readConfig(config, "GM",                GM,    Config::DEFAULT,  STRING_DEFAULT_GM, "Geocentric gravitational constant");
      if(!isCreateSchema(config))
      {
        Kepler kepler(Time(), Omega, i, omega, a, e, M, GM);
        kepler.orbit(Time(), startEpoch.position, startEpoch.velocity);
      }
    }
    if(readConfigChoiceElement(config, "positionAndVelocity", choice, ""))
    {
      readConfig(config, "position0x", startEpoch.position.x(), Config::MUSTSET, "", "[m] in CRF");
      readConfig(config, "position0y", startEpoch.position.y(), Config::MUSTSET, "", "[m] in CRF");
      readConfig(config, "position0z", startEpoch.position.z(), Config::MUSTSET, "", "[m] in CRF");
      readConfig(config, "velocity0x", startEpoch.velocity.x(), Config::MUSTSET, "", "[m/s]");
      readConfig(config, "velocity0y", startEpoch.velocity.y(), Config::MUSTSET, "", "[m/s]");
      readConfig(config, "velocity0z", startEpoch.velocity.z(), Config::MUSTSET, "", "[m/s]");
    }
    if(readConfigChoiceElement(config, "file", choice, ""))
    {
      readConfig(config, "inputfileOrbit", fileNameStart, Config::MUSTSET,  "",     "only epoch at timeStart is used");
      readConfig(config, "margin",         margin,        Config::DEFAULT,  "1e-5", "[seconds] used when finding initial epoch in orbitFile");
    }
    endChoice(config);

    readConfig(config, "propagator",     orbitPropagator,   Config::MUSTSET,  "polynomial", "orbit propagation method.");
    readConfig(config, "earthRotation",  earthRotation,     Config::MUSTSET,  "",    "");
    readConfig(config, "ephemerides",    ephemerides,       Config::OPTIONAL, "jpl", "");
    readConfig(config, "forces",         forces,            Config::MUSTSET,  "",    "considered in orbit propagation.");
    readConfig(config, "reverse",        reverse,           Config::DEFAULT,  "0",   "start integration at last epoch in timeSeries, going backward in time.");
    if(isCreateSchema(config)) return;

    // Setup
    // -----
    SatelliteModelPtr satellite;
    if(!fileNameSatellite.empty())
      readFileSatelliteModel(fileNameSatellite, satellite);

    std::vector<Time> times = timeSeries->times();
    Time timeStart = times.front();
    Time sampling  = medianSampling(times);
    if(!isRegular(times))
      throw(Exception("Time intervals must be regularly spaced."));

    // Need to flip?
    if(reverse)
    {
      timeStart = times.back();
      sampling  = -sampling;
    }

    // read startEpoch from file
    // -------------------------
    if(!fileNameStart.empty())
    {
      OrbitArc orbit = InstrumentFile::read(fileNameStart);
      for(UInt i=0; i<orbit.size(); i++)
        if(std::fabs((orbit.at(i).time-timeStart).seconds())<margin)
          startEpoch = orbit.at(i);
      if(std::fabs((startEpoch.time-timeStart).seconds())>margin)
        throw (Exception("Requested start time "+timeStart.dateTimeStr()+" not present in file <"+fileNameStart.str()+">"));
    }

    // Integrate
    // ---------
    logStatus<<"integrating orbit"<<Log::endl;
    startEpoch.time = timeStart;
    OrbitArc orbit = orbitPropagator->integrateArc(startEpoch, sampling, times.size(), forces, satellite, earthRotation, ephemerides);
    if(reverse)
      orbit = OrbitPropagator::flip(orbit);

    // Save
    // ----
    logStatus<<"write orbit data to file <"<<fileNameOrbit<<">"<<Log::endl;
    InstrumentFile::write(fileNameOrbit, orbit);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
