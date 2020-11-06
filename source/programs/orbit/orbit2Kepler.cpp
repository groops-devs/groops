/***********************************************/
/**
* @file orbit2Kepler.cpp
*
* @brief Keplerian elements from orbit position and velocity at each epoch.
*
* @author Torsten Mayer-Guerr
* @date 2010-07-20
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program computes the osculating Keplerian elements from position and velocity
of a given \configFile{inputfileOrbit}{instrument}.
The \configFile{inputfileOrbit}{instrument} must contain positions and velocities (see \program{OrbitAddVelocityAndAcceleration}).

The \config{outputfileKepler} is an \file{instrument file}{instrument} (MISCVALUES)
with the Keplerian elements at each epoch in the following order
\begin{itemize}
\item Ascending Node $\Omega$ [degree]
\item Inclination $i$ [degree]
\item Argument of perigee $\omega$ [degree]
\item major axis $a$ [m]
\item eccentricity $e$
\item mean anomaly $M$ [degree]
\item transit time of perigee $\tau$ [mjd]
\end{itemize}
)";

/***********************************************/

#include "programs/program.h"
#include "base/kepler.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Keplerian elements from orbit position and velocity at each epoch.
* @ingroup programsGroup */
class Orbit2Kepler
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(Orbit2Kepler, PARALLEL, "keplerian elements from orbit position and velocity at each epoch", Orbit, Instrument, TimeSeries)

/***********************************************/

void Orbit2Kepler::run(Config &config)
{
  try
  {
    FileName              fileNameOut, fileNameOrbit;
    std::vector<FileName> fileNamesInstrument;
    Double                GM;

    readConfig(config, "outputfileKepler",    fileNameOut,         Config::MUSTSET,  "", "each epoch: Omega, i, omega [degree], a [m], e, M [degree], tau [mjd]");
    readConfig(config, "inputfileOrbit",      fileNameOrbit,       Config::MUSTSET,  "", "position and velocity at each epoch define the kepler orbit");
    readConfig(config, "inputfileInstrument", fileNamesInstrument, Config::OPTIONAL, "", "data is appended");
    readConfig(config, "GM",                  GM,                  Config::DEFAULT,  STRING_DEFAULT_GM, "Geocentric gravitational constant");
    if(isCreateSchema(config)) return;

    // =======================

    logStatus<<"calculate Keplerian elements"<<Log::endl;
    InstrumentFile orbitFile(fileNameOrbit);
    UInt dataCount = 8;
    std::vector<InstrumentFilePtr> instrumentFile;
    for(auto &fileName : fileNamesInstrument)
    {
      instrumentFile.push_back(InstrumentFile::newFile(fileName));
      InstrumentFile::checkArcCount({orbitFile, *instrumentFile.back()});
      dataCount += instrumentFile.back()->dataCount(TRUE/*mustDefined*/);
    }

    std::vector<Arc> arcList(orbitFile.arcCount());
    Parallel::forEach(arcList, [&] (UInt arcNo)
    {
      const OrbitArc orbit = orbitFile.readArc(arcNo);
      Matrix A(orbit.size(), dataCount);
      for(UInt i=0; i<orbit.size(); i++)
      {
        if(orbit.at(i).velocity.r()==0)
          throw(Exception("no velocity given"));

        const Kepler kepler(orbit.at(i).time, orbit.at(i).position, orbit.at(i).velocity, GM);
        const Double M = (orbit.at(i).time-kepler.tau).seconds() * std::sqrt(GM/std::pow(kepler.a,3));

        A(i, 1) = std::fmod(kepler.Omega+2*PI, 2*PI) * DEG2RAD;
        A(i, 2) = kepler.i * DEG2RAD;
        A(i, 3) = std::fmod(kepler.omega+2*PI, 2*PI) * DEG2RAD;
        A(i, 4) = kepler.a;
        A(i, 5) = kepler.e;
        A(i, 6) = std::fmod(M+2*PI, 2*PI) * DEG2RAD;
        A(i, 7) = kepler.tau.mjd();
      }

      UInt idx = 8;
      for(auto &file: instrumentFile)
      {
        Arc arc = file->readArc(arcNo);
        Arc::checkSynchronized({orbit, arc});
        Matrix B = arc.matrix();
        copy(B.column(1, B.columns()-1), A.column(idx, B.columns()-1));
        idx += B.columns()-1;
      }

      return Arc(orbit.times(), A);
    });

    // write results
    // -------------
    if(Parallel::isMaster())
    {
      logStatus<<"write time series of Keplerian elements to file <"<<fileNameOut<<">"<<Log::endl;
      InstrumentFile::write(fileNameOut, arcList);
      Arc::printStatistics(arcList);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
