/***********************************************/
/**
* @file orbit2EclipseFactor.cpp
*
* @brief Create instrument file containing eclipse factors.
*
* @author Beate Klinger
* @date 2016-06-10
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program generates an instrument file containing the eclipse factor for a given set of orbit.
The data of \configFile{inputfileInstrument}{instrument} are appended as values to each point.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "classes/ephemerides/ephemerides.h"
#include "classes/eclipse/eclipse.h"

/***** CLASS ***********************************/

/** @brief Create instrument file containing eclipse factors.
* @ingroup programsGroup */
class Orbit2EclipseFactor
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(Orbit2EclipseFactor, PARALLEL, "Create instrument file containing eclipse factors.", Orbit, Instrument)
GROOPS_RENAMED_PROGRAM(InstrumentOrbit2EclipseFactor, Orbit2EclipseFactor, date2time(2020, 05, 25))

/***********************************************/

void Orbit2EclipseFactor::run(Config &config)
{
  try
  {
    FileName              fileNameEclipse, fileNameOrbit;
    std::vector<FileName> fileNamesInstrument;
    EphemeridesPtr        ephemerides;
    EclipsePtr            eclipse;

    readConfig(config, "outputfileEclipseFactor", fileNameEclipse,     Config::MUSTSET,  "", "");
    readConfig(config, "inputfileOrbit",          fileNameOrbit,       Config::MUSTSET,  "", "");
    readConfig(config, "inputfileInstrument",     fileNamesInstrument, Config::OPTIONAL, "", "data are appended");
    readConfig(config, "ephemerides",             ephemerides,         Config::MUSTSET,  "", "");
    readConfig(config, "eclipse",                 eclipse,             Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    // =======================

    logStatus<<"computing eclipse factor"<<Log::endl;
    InstrumentFile orbitFile(fileNameOrbit);
    UInt dataCount = 2;
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
        A(i, 1) = eclipse->factor(orbit.at(i).time, orbit.at(i).position, ephemerides);

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
      logStatus<<"write eclipse factor to file <"<<fileNameEclipse<<">"<<Log::endl;
      InstrumentFile::write(fileNameEclipse, arcList);
      Arc::printStatistics(arcList);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
