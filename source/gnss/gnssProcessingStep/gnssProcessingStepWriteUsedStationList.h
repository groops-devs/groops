/***********************************************/
/**
* @file gnssProcessingStepWriteUsedStationList.h
*
* @brief GNSS processing step: WriteUsedStationList.
*
* @author Torsten Mayer-Guerr
* @date 2021-09-05
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPROCESSINGSTEPWRITEUSEDSTATIONLIST__
#define __GROOPS_GNSSPROCESSINGSTEPWRITEUSEDSTATIONLIST__

// Latex documentation
#ifdef DOCSTRING_GnssProcessingStep
static const char *docstringGnssProcessingStepWriteUsedStationList = R"(
\subsection{WriteUsedStationList}\label{gnssProcessingStepType:writeUsedStationList}
Writes a \file{list}{stringList} of receivers (stations) which are used in the last step and
selected by \configClass{selectReceivers}{gnssTransceiverSelectorType}.
)";
#endif

/***********************************************/

#include "config/config.h"
#include "files/fileStringTable.h"
#include "gnss/gnssTransceiverSelector/gnssTransceiverSelector.h"
#include "gnss/gnssProcessingStep/gnssProcessingStep.h"

/***** CLASS ***********************************/

/** @brief GNSS processing step: WriteUsedStationList.
* @ingroup gnssProcessingStepGroup
* @see GnssProcessingStep */
class GnssProcessingStepWriteUsedStationList : public GnssProcessingStepBase
{
  GnssTransceiverSelectorPtr selectReceivers;
  FileName                   fileNameUsedStationList;

public:
  GnssProcessingStepWriteUsedStationList(Config &config);
  void process(GnssProcessingStep::State &state) override;
};

/***********************************************/

inline GnssProcessingStepWriteUsedStationList::GnssProcessingStepWriteUsedStationList(Config &config)
{
  try
  {
    readConfig(config, "selectReceivers",           selectReceivers,         Config::MUSTSET, "", "subset of used stations");
    readConfig(config, "outputfileUsedStationList", fileNameUsedStationList, Config::MUSTSET, "", "ascii file with names of used stations");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssProcessingStepWriteUsedStationList::process(GnssProcessingStep::State &state)
{
  try
  {
    if(!Parallel::isMaster(state.normalEquationInfo.comm))
      return;
    auto selectedReceivers = selectReceivers->select(state.gnss->receivers);
    logStatus<<"write used station list to file <"<<fileNameUsedStationList<<">"<<Log::endl;
    std::vector<std::string> usedStationList;
    for(auto recv : state.gnss->receivers)
      if(selectedReceivers.at(recv->idRecv()) && state.normalEquationInfo.estimateReceiver.at(recv->idRecv()))
        usedStationList.push_back(recv->name());
    writeFileStringList(fileNameUsedStationList, usedStationList);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
