/***********************************************/
/**
* @file instrumentApplyTimeOffset.cpp
*
* @brief Apply a time offset to an instrument file.
*
* @author Andreas Kvas
* @author Torsten Mayer-Guerr
* @date 2022-04-01
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program applies a \configFile{inputfileTimeOffset}{instrument} (MISCVALUE)
to an \configFile{inputfileInstrument}{instrument}.
The time offsets in seconds are multiplicated with a \config{factor}.
The instrument files must be synchronized (see \program{InstrumentSynchronize}).
)";

/***********************************************/

#include "programs/program.h"
#include "base/polynomial.h"
#include "files/fileInstrument.h"
#include "files/fileMatrix.h"

/***** CLASS ***********************************/

/** @brief Apply a time offset to an instrument file.
* @ingroup programsGroup */
class InstrumentApplyTimeOffset
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(InstrumentApplyTimeOffset, SINGLEPROCESS, "Apply a time offset to an instrument file", Instrument)

/***********************************************/

void InstrumentApplyTimeOffset::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName fileNameOut, fileNameIn, fileNameTimeOffset;
    Double   factor;

    readConfig(config, "outputfileInstrument", fileNameOut,        Config::MUSTSET, "",    "");
    readConfig(config, "inputfileInstrument",  fileNameIn,         Config::MUSTSET, "",    "");
    readConfig(config, "inputfileTimeOffset",  fileNameTimeOffset, Config::MUSTSET, "",    "MISCVALUE with time offset in seconds");
    readConfig(config, "factor",               factor,             Config::DEFAULT, "1.0", "applied to time offset");
    if(isCreateSchema(config)) return;

    logStatus<<"read instrument data <"<<fileNameIn<<">"<<Log::endl;
    InstrumentFile instrumentFile(fileNameIn);
    InstrumentFile timeOffsetFile(fileNameTimeOffset);
    InstrumentFile::checkArcCount({instrumentFile, timeOffsetFile});

    std::vector<Arc> arcList(instrumentFile.arcCount());
    Parallel::forEach(arcList, [&](UInt arcNo)
    {
      Arc          arc           = instrumentFile.readArc(arcNo);
      MiscValueArc arcTimeOffset = timeOffsetFile.readArc(arcNo);
      Arc::checkSynchronized({arc, arcTimeOffset});

      for(UInt i=0; i<arc.size(); i++)
        arc.at(i).time += factor * seconds2time(arcTimeOffset.at(i).value);
      return arc;
    }, comm);

    if(Parallel::isMaster(comm))
    {
      logStatus<<"write instrument data to file <"<<fileNameOut<<">"<<Log::endl;
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
