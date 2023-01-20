/***********************************************/
/**
* @file gnssProcessingStepWriteUsedTransmitterList.h
*
* @brief GNSS processing step: WriteUsedTransmitterList.
*
* @author Torsten Mayer-Guerr
* @date 2021-09-05
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPROCESSINGSTEPWRITEUSEDTRANSMITTERLIST__
#define __GROOPS_GNSSPROCESSINGSTEPWRITEUSEDTRANSMITTERLIST__

// Latex documentation
#ifdef DOCSTRING_GnssProcessingStep
static const char *docstringGnssProcessingStepWriteUsedTransmitterList = R"(
\subsection{WriteUsedTransmitterList}\label{gnssProcessingStepType:writeUsedTransmitterList}
Writes a \file{list}{stringList} of transmitters which are used in the last step and
selected by \configClass{selectTransmitters}{platformSelectorType}.
)";
#endif

/***********************************************/

#include "config/config.h"
#include "files/fileStringTable.h"
#include "classes/platformSelector/platformSelector.h"
#include "gnss/gnssProcessingStep/gnssProcessingStep.h"

/***** CLASS ***********************************/

/** @brief GNSS processing step: WriteUsedTransmitterList.
* @ingroup gnssProcessingStepGroup
* @see GnssProcessingStep */
class GnssProcessingStepWriteUsedTransmitterList : public GnssProcessingStepBase
{
  PlatformSelectorPtr selectTransmitters;
  FileName            fileNameOutTransmitterList;

public:
  GnssProcessingStepWriteUsedTransmitterList(Config &config);
  void process(GnssProcessingStep::State &state) override;
};

/***********************************************/

inline GnssProcessingStepWriteUsedTransmitterList::GnssProcessingStepWriteUsedTransmitterList(Config &config)
{
  try
  {
    readConfig(config, "selectTransmitters",            selectTransmitters,         Config::MUSTSET, "", "subset of used transmitters");
    readConfig(config, "outputfileUsedTransmitterList", fileNameOutTransmitterList, Config::MUSTSET, "", "ascii file with PRNs");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void GnssProcessingStepWriteUsedTransmitterList::process(GnssProcessingStep::State &state)
{
  try
  {
    if(!Parallel::isMaster(state.normalEquationInfo.comm))
      return;
    auto selectedTransmitters = state.gnss->selectTransmitters(selectTransmitters);
    logStatus<<"write used transmitter list to file <"<<fileNameOutTransmitterList<<">"<<Log::endl;
    std::vector<std::string> usedTransmitterList;
    for(auto trans : state.gnss->transmitters)
      if(selectedTransmitters.at(trans->idTrans()) && trans->useable())
        usedTransmitterList.push_back(trans->name());
    writeFileStringList(fileNameOutTransmitterList, usedTransmitterList);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
