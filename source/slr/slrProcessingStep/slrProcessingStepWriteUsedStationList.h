/***********************************************/
/**
* @file slrProcessingStepWriteUsedStationList.h
*
* @brief SLR processing step: WriteUsedStationList.
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#ifndef __GROOPS_SLRPROCESSINGSTEPWRITEUSEDSTATIONLIST__
#define __GROOPS_SLRPROCESSINGSTEPWRITEUSEDSTATIONLIST__

// Latex documentation
#ifdef DOCSTRING_SlrProcessingStep
static const char *docstringSlrProcessingStepWriteUsedStationList = R"(
\subsection{WriteUsedStationList}\label{slrProcessingStepType:writeUsedStationList}
Writes a \file{list}{stringList} of stations (stations) which are used in the last step and
selected by \configClass{selectStations}{platformSelectorType}.
)";
#endif

/***********************************************/

#include "config/config.h"
#include "files/fileStringTable.h"
#include "classes/platformSelector/platformSelector.h"
#include "slr/slrProcessingStep/slrProcessingStep.h"

/***** CLASS ***********************************/

/** @brief SLR processing step: WriteUsedStationList.
* @ingroup slrProcessingStepGroup
* @see SlrProcessingStep */
class SlrProcessingStepWriteUsedStationList : public SlrProcessingStepBase
{
  PlatformSelectorPtr selectorStations;
  FileName            fileNameUsedStationList;

public:
  SlrProcessingStepWriteUsedStationList(Config &config);
  void process(SlrProcessingStep::State &state) override;
};

/***********************************************/

inline SlrProcessingStepWriteUsedStationList::SlrProcessingStepWriteUsedStationList(Config &config)
{
  try
  {
    readConfig(config, "selectStations",            selectorStations,        Config::MUSTSET, "", "subset of used stations");
    readConfig(config, "outputfileUsedStationList", fileNameUsedStationList, Config::MUSTSET, "output/stationList.txt", "ascii file with names of used stations");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void SlrProcessingStepWriteUsedStationList::process(SlrProcessingStep::State &state)
{
  try
  {
    auto selectedStations = state.slr->selectStations(selectorStations);
    logStatus<<"write used station list to file <"<<fileNameUsedStationList<<">"<<Log::endl;
    std::vector<std::string> usedStationList;
    for(auto stat : state.slr->stations)
      if(selectedStations.at(stat->idStat()) && state.normalEquationInfo.estimateStation.at(stat->idStat()) && stat->useable())
        usedStationList.push_back(stat->name());
    writeFileStringList(fileNameUsedStationList, usedStationList);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
