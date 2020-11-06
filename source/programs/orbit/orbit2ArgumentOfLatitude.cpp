/***********************************************/
/**
* @file orbit2ArgumentOfLatitude.cpp
*
* @brief Computes argument of latitude.
*
* @author Torsten Mayer-Guerr
* @date 2016-07-26
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program computes the argument of latitude of an \file{orbit}{instrument}
and writes it as \file{instrument file}{instrument} (MISCVALUES).
The data of \configFile{inputfileInstrument}{instrument} are appended as values to each point.

\fig{!hb}{0.8}{instrumentOrbit2ArgumentOfLatitude}{fig:instrumentOrbit2ArgumentOfLatitude}{Derivation filtered GRACE range-rate residuals.}
)";

/***********************************************/

#include "programs/program.h"
#include "base/planets.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Argument of latitude.
* @ingroup programsGroup */
class Orbit2ArgumentOfLatitude
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(Orbit2ArgumentOfLatitude, PARALLEL, "Argument of latitude.", Orbit, Instrument)
GROOPS_RENAMED_PROGRAM(InstrumentOrbit2ArgumentOfLatitude, Orbit2ArgumentOfLatitude, date2time(2020, 05, 25))

/***********************************************/

void Orbit2ArgumentOfLatitude::run(Config &config)
{
  try
  {
    FileName fileNameOut, fileNameOrbit;
    std::vector<FileName> fileNamesInstrument;

    readConfig(config, "outputfileArgOfLatitude", fileNameOut,         Config::MUSTSET,  "", "instrument file (MISCVALUES)");
    readConfig(config, "inputfileOrbit",          fileNameOrbit,       Config::MUSTSET,  "", "");
    readConfig(config, "inputfileInstrument",     fileNamesInstrument, Config::OPTIONAL, "", "data are appended");
    if(isCreateSchema(config)) return;

    // =======================

    logStatus<<"computing argument of latitude"<<Log::endl;
    InstrumentFile orbitFile(fileNameOrbit);
    UInt dataCount = 2;  // time, arg
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
        const Vector3d z = normalize(crossProduct(orbit.at(i).position, orbit.at(i).velocity));
        const Vector3d x = normalize(crossProduct(Vector3d(0,0,1), z));
        const Vector3d y = crossProduct(z, x);
        A(i, 1) = atan2(inner(orbit.at(i).position,y), inner(orbit.at(i).position, x)); // Argument of latitude of satellite
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
      logStatus<<"write argument of latitude to file <"<<fileNameOut<<">"<<Log::endl;
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
