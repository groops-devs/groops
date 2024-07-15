/***********************************************/
/**
* @file slrProcessingStepSelectStations.h
*
* @brief SLR processing step: SelectStations.
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#ifndef __GROOPS_SLRPROCESSINGSTEPSELECTSTATION__
#define __GROOPS_SLRPROCESSINGSTEPSELECTSTATION__

// Latex documentation
#ifdef DOCSTRING_SlrProcessingStep
static const char *docstringSlrProcessingStepSelectStations = R"(
\subsection{SelectStations}\label{slrProcessingStepType:selectStations}
This step can be used to process only a subset of stations in subsequent processing steps.
)";
#endif

/***********************************************/

#include "config/config.h"
#include "classes/platformSelector/platformSelector.h"
#include "slr/slrProcessingStep/slrProcessingStep.h"

/***** CLASS ***********************************/

/** @brief SLR processing step: SelectStations.
* @ingroup slrProcessingStepGroup
* @see SlrProcessingStep */
class SlrProcessingStepSelectStations : public SlrProcessingStepBase
{
  PlatformSelectorPtr selectorStations;

public:
  SlrProcessingStepSelectStations(Config &config);
  void process(SlrProcessingStep::State &state) override;
  Bool expectInitializedParameters() const override {return FALSE;}
};

/***********************************************/

inline SlrProcessingStepSelectStations::SlrProcessingStepSelectStations(Config &config)
{
  try
  {
    readConfig(config, "selectStations", selectorStations, Config::MUSTSET,  "", "");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void SlrProcessingStepSelectStations::process(SlrProcessingStep::State &state)
{
  try
  {
    logStatus<<"=== select stations ========================================"<<Log::endl;
    state.normalEquationInfo.estimateStation = state.slr->selectStations(selectorStations);
    for(auto stat : state.slr->stations)
      if(!stat->useable())
        state.normalEquationInfo.estimateStation.at(stat->idStat()) = FALSE;
    logInfo<<"  "<<std::count(state.normalEquationInfo.estimateStation.begin(), state.normalEquationInfo.estimateStation.end(), TRUE)<<" stations selected"<<Log::endl;
    state.changedNormalEquationInfo = TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
