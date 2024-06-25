/***********************************************/
/**
* @file slrProcessingStepWriteUsedSatelliteList.h
*
* @brief SLR processing step: WriteUsedSatelliteList.
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#ifndef __GROOPS_SLRPROCESSINGSTEPWRITEUSEDSATELLITELIST__
#define __GROOPS_SLRPROCESSINGSTEPWRITEUSEDSATELLITELIST__

// Latex documentation
#ifdef DOCSTRING_SlrProcessingStep
static const char *docstringSlrProcessingStepWriteUsedSatelliteList = R"(
\subsection{WriteUsedSatelliteList}\label{slrProcessingStepType:writeUsedSatelliteList}
Writes a \file{list}{stringList} of satellites which are used in the last step and
selected by \configClass{selectSatellites}{platformSelectorType}.
)";
#endif

/***********************************************/

#include "config/config.h"
#include "files/fileStringTable.h"
#include "classes/platformSelector/platformSelector.h"
#include "slr/slrProcessingStep/slrProcessingStep.h"

/***** CLASS ***********************************/

/** @brief SLR processing step: WriteUsedSatelliteList.
* @ingroup slrProcessingStepGroup
* @see SlrProcessingStep */
class SlrProcessingStepWriteUsedSatelliteList : public SlrProcessingStepBase
{
  PlatformSelectorPtr selectorSatellites;
  FileName            fileNameOutSatelliteList;

public:
  SlrProcessingStepWriteUsedSatelliteList(Config &config);
  void process(SlrProcessingStep::State &state) override;
};

/***********************************************/

inline SlrProcessingStepWriteUsedSatelliteList::SlrProcessingStepWriteUsedSatelliteList(Config &config)
{
  try
  {
    readConfig(config, "selectSatellites",            selectorSatellites,       Config::MUSTSET, "", "subset of used satellites");
    readConfig(config, "outputfileUsedSatelliteList", fileNameOutSatelliteList, Config::MUSTSET, "output/satelliteList.txt", "ascii file with names");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void SlrProcessingStepWriteUsedSatelliteList::process(SlrProcessingStep::State &state)
{
  try
  {
    auto selectedSatellites = state.slr->selectSatellites(selectorSatellites);
    logStatus<<"write used satellite list to file <"<<fileNameOutSatelliteList<<">"<<Log::endl;
    std::vector<std::string> usedSatelliteList;
    for(auto sat : state.slr->satellites)
      if(selectedSatellites.at(sat->idSat()) && state.normalEquationInfo.estimateSatellite.at(sat->idSat()) && sat->useable())
        usedSatelliteList.push_back(sat->name());
    writeFileStringList(fileNameOutSatelliteList, usedSatelliteList);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
