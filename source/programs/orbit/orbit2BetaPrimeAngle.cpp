/***********************************************/
/**
* @file orbit2BetaPrimeAngle.cpp
*
* @brief Computes beta prime angle.
*
* @author Beate Klinger
* @author Norbert Zehentner
* @date 2016-08-01
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program computes the beta prime angle (between the orbital plane and earth-sun direction)
and writes it as MISCVALUES \file{instrument file}{instrument}. The angle is calculated w.r.t the sun (per default),
but can be changed.
The data of \configFile{inputfileInstrument}{instrument} are appended as values to each point.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "classes/ephemerides/ephemerides.h"

/***** CLASS ***********************************/

/** @brief Beta prime angle.
* @ingroup programsGroup */
class Orbit2BetaPrimeAngle
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(Orbit2BetaPrimeAngle, PARALLEL, "Beta prime angle", Orbit, Instrument)
GROOPS_RENAMED_PROGRAM(InstrumentOrbit2BetaPrimeAngle, Orbit2BetaPrimeAngle, date2time(2020, 05, 25))

/***********************************************/

void Orbit2BetaPrimeAngle::run(Config &config)
{
  try
  {
    FileName              fileNameOut, fileNameOrbit;
    std::vector<FileName> fileNamesInstrument;
    EphemeridesPtr        ephemerides;
    Ephemerides::Planet   planet;
    std::string           choice;

    readConfig(config, "outputfileBetaAngle", fileNameOut,         Config::MUSTSET,  "",    "instrument file (MISCVALUES)");
    readConfig(config, "inputfileOrbit",      fileNameOrbit,       Config::MUSTSET,  "",    "");
    readConfig(config, "inputfileInstrument", fileNamesInstrument, Config::OPTIONAL, "",    "data are appended");
    readConfig(config, "ephemerides",         ephemerides,         Config::MUSTSET,  "",    "");
    readConfig(config, "planet",              planet,              Config::DEFAULT,  "sun", "");
    if(isCreateSchema(config)) return;

    // =======================

    logStatus<<"computing beta prime angle"<<Log::endl;
    InstrumentFile orbitFile(fileNameOrbit);
    UInt dataCount = 2;  // time, beta
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
        const Vector3d planetPos = normalize(ephemerides->position(orbit.at(i).time, planet)); // direction from central-body to planet (earth-sun direction)
        A(i, 1) = std::acos(inner(normalize(crossProduct(normalize(orbit.at(i).velocity), normalize(orbit.at(i).position))), planetPos)) - PI/2;
      }

      UInt idx = 2;
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
      logStatus<<"write beta prime angle to file <"<<fileNameOut<<">"<<Log::endl;
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
