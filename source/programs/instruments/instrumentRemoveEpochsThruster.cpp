/***********************************************/
/**
* @file instrumentRemoveEpochsThruster.cpp
*
* @brief Remove epochs from instrument data.
*
* @author Torsten Mayer-Guerr
* @date 2022-07-30
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program remove epochs from an \file{instrument file}{instrument}.
The epochs are defined by a \file{thruster file}{instrument}
plus a defined margin before and after the thruster firings.
The arcs of the instrument file are concatenated to one arc.
The removed epochs can be saved in a separate instrument file.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "classes/timeSeries/timeSeries.h"

/***** CLASS ***********************************/

/** @brief Remove epochs from instrument data.
* @ingroup programsGroup */
class InstrumentRemoveEpochsThruster
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(InstrumentRemoveEpochsThruster, SINGLEPROCESS, "Remove epochs from instrument data", Instrument)

/***********************************************/

void InstrumentRemoveEpochsThruster::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameOut, fileNameOutRemoved, fileNameIn, fileNameThruster;
    Double   marginBefore, marginAfter;

    readConfig(config, "outputfileInstrument",               fileNameOut,        Config::OPTIONAL, "",    "all epochs are concatenated in one arc");
    readConfig(config, "outputfileInstrumentRemovedEpochs",  fileNameOutRemoved, Config::OPTIONAL, "",    "all epochs are concatenated in one arc");
    readConfig(config, "inputfileInstrument",                fileNameIn,         Config::MUSTSET,  "",    "");
    readConfig(config, "inputfileThruster",                  fileNameThruster,   Config::MUSTSET,  "",    "THRUSTER");
    readConfig(config, "marginBefore",                       marginBefore,       Config::DEFAULT,  "0.1", "margin before start of firing [seconds]");
    readConfig(config, "marginAfter",                        marginAfter,        Config::DEFAULT,  "0.7", "margin after end of firing [seconds]");
    if(isCreateSchema(config)) return;

    // thruster time series
    // --------------------
    logStatus<<"read thruster data <"<<fileNameThruster<<">"<<Log::endl;
    ThrusterArc thruster = InstrumentFile::read(fileNameThruster);
    std::vector<Time> timesStart = thruster.times();
    std::vector<Time> timesEnd   = timesStart;
    for(UInt i=0; i<thruster.size(); i++)
      timesStart.at(i) -= seconds2time(marginBefore);
    for(UInt i=0; i<thruster.size(); i++)
      timesEnd.at(i) += seconds2time(1e-3*max(thruster.at(i).data()) + marginAfter);

    logStatus<<"read instrument data <"<<fileNameIn<<">"<<Log::endl;
    Arc arc = InstrumentFile::read(fileNameIn);

    logStatus<<"remove thruster epochs from instrument data"<<Log::endl;
    Arc arcNew, arcRemoved;
    UInt idxThruster = 0;
    Single::forEach(arc.size(), [&](UInt i)
    {
      while((idxThruster < timesEnd.size()) && ((arc.at(i).time-timesEnd.at(idxThruster)).seconds() >= 0)) // find next thruster epoch
        idxThruster++;
      if((idxThruster >= timesEnd.size()) || ((arc.at(i).time-timesStart.at(idxThruster)).seconds() < 0))
        arcNew.push_back(arc.at(i));
      else
        arcRemoved.push_back(arc.at(i));
    });

    // save instrument data
    // --------------------
    if(!fileNameOut.empty())
    {
      logStatus<<"write instrument data to file <"<<fileNameOut<<">"<<Log::endl;
      InstrumentFile::write(fileNameOut, arcNew);
      Arc::printStatistics(arcNew);
    }

    if(!fileNameOutRemoved.empty())
    {
      logStatus<<"write removed epochs to file <"<<fileNameOutRemoved<<">"<<Log::endl;
      InstrumentFile::write(fileNameOutRemoved, arcRemoved);
      Arc::printStatistics(arcRemoved);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
