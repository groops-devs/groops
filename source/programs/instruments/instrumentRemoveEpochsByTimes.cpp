/***********************************************/
/**
* @file instrumentRemoveEpochsByTimes.cpp
*
* @brief Remove epochs from instrument data.
*
* @author Beate Klinger
* @date 2015-02-03
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program compares an \file{instrument file}{instrument} with a
\configClass{time series}{timeSeriesType}.
Epochs contained within the time series (including a defined margin)
are removed from the instrument file. The margin is added on
both sides of the epochs. The arcs of the instrument file are
concatenated to one arc. The removed epochs can be saved
in a separate instrument file.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "classes/timeSeries/timeSeries.h"

/***** CLASS ***********************************/

/** @brief Remove epochs from instrument data.
* @ingroup programsGroup */
class InstrumentRemoveEpochsByTimes
{
  InstrumentFile timeFile;

public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(InstrumentRemoveEpochsByTimes, SINGLEPROCESS, "Remove epochs from instrument data", Instrument)

/***********************************************/

void InstrumentRemoveEpochsByTimes::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName      outName, outNameRemoved, inName;
    Double        margin;
    TimeSeriesPtr timeSeries;

    readConfig(config, "outputfileInstrument",               outName,        Config::OPTIONAL, "", "all epochs are concatenated in one arc");
    readConfig(config, "outputfileInstrumentRemovedEpochs",  outNameRemoved, Config::OPTIONAL, "", "all epochs are concatenated in one arc");
    readConfig(config, "inputfileInstrument",                inName,         Config::MUSTSET,  "", "");
    readConfig(config, "timePoints",                         timeSeries,     Config::MUSTSET,  "", "");
    readConfig(config, "margin",                             margin,         Config::DEFAULT,  "1e-5", "margin size (on both sides) [seconds]");
    if(isCreateSchema(config)) return;

    // ======================================================

    // read time series
    // ----------------
    logStatus<<"read time series"<<Log::endl;
    std::vector<Time> times = timeSeries->times();
    logStatus<<"  epochs:  "<<times.size()<<Log::endl;

    // quick test
    // ----------
    if((times.size()==0) && (inName.str() == outName.str()))
    {
      logStatus<<"No epochs to remove."<<Log::endl;
      return;
    }

    // read instrument data
    // --------------------
    logStatus<<"read instrument data <"<<inName<<">"<<Log::endl;
    Arc arc = InstrumentFile::read(inName);

    // ======================================================

    // remove epochs within buffer
    // ---------------------------
    Arc arcNew, arcRemoved;
    if(times.size())
    {
      logStatus<<"remove epochs (+/- "<<margin<<" sec) from instrument data"<<Log::endl;
      UInt idxTime = 0;
      Single::forEach(arc.size(), [&](UInt i)
      {
        while((idxTime < times.size()) && ((arc.at(i).time-times.at(idxTime)).seconds() > margin))
          idxTime++;

        if((idxTime < times.size()) && ((times.at(idxTime)-arc.at(i).time).seconds() <= margin))
          arcRemoved.push_back(arc.at(i));
        else
          arcNew.push_back(arc.at(i));
      });
      logInfo<<"  "<<arcRemoved.size()<<" epochs removed"<<Log::endl;
    }
    else
    {
      logStatus<<"No epochs to remove."<<Log::endl;
      arcNew = arc;
    }

    // ======================================================

    // save instrument data
    // --------------------
    if(!outName.empty())
    {
      logStatus<<"write instrument data to file <"<<outName<<">"<<Log::endl;
      InstrumentFile::write(outName, arcNew);
      Arc::printStatistics(arcNew);
    }

    if(!outNameRemoved.empty())
    {
      logStatus<<"write removed epochs to file <"<<outNameRemoved<<">"<<Log::endl;
      InstrumentFile::write(outNameRemoved, arcRemoved);
      Arc::printStatistics(arcRemoved);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
