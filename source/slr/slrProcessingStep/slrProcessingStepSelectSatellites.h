/***********************************************/
/**
* @file slrProcessingStepSelectSatellites.h
*
* @brief SLR processing step: SelectSatellites.
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#ifndef __GROOPS_SLRPROCESSINGSTEPSELECTSATELLITES__
#define __GROOPS_SLRPROCESSINGSTEPSELECTSATELLITES__

// Latex documentation
#ifdef DOCSTRING_SlrProcessingStep
static const char *docstringSlrProcessingStepSelectSatellites = R"(
\subsection{SelectSatellites}\label{slrProcessingStepType:selectSatellites}
This step can be used to process only a subset of satellites in subsequent processing steps.
)";
#endif

/***********************************************/

#include "config/config.h"
#include "classes/platformSelector/platformSelector.h"
#include "slr/slrProcessingStep/slrProcessingStep.h"

/***** CLASS ***********************************/

/** @brief SLR processing step: SelectSatellites.
* @ingroup slrProcessingStepGroup
* @see SlrProcessingStep */
class SlrProcessingStepSelectSatellites : public SlrProcessingStepBase
{
  PlatformSelectorPtr selectorSatellites;

public:
  SlrProcessingStepSelectSatellites(Config &config);
  void process(SlrProcessingStep::State &state) override;
  Bool expectInitializedParameters() const override {return FALSE;}
};

/***********************************************/

inline SlrProcessingStepSelectSatellites::SlrProcessingStepSelectSatellites(Config &config)
{
  try
  {
    readConfig(config, "selectSatellites", selectorSatellites, Config::MUSTSET,  "", "");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void SlrProcessingStepSelectSatellites::process(SlrProcessingStep::State &state)
{
  try
  {
    logStatus<<"=== select satellites ========================================"<<Log::endl;
    state.normalEquationInfo.estimateSatellite = state.slr->selectSatellites(selectorSatellites);
    for(auto sat : state.slr->satellites)
      if(!sat->useable())
        state.normalEquationInfo.estimateSatellite.at(sat->idSat()) = FALSE;
    logInfo<<"  "<<std::count(state.normalEquationInfo.estimateSatellite.begin(), state.normalEquationInfo.estimateSatellite.end(), TRUE)<<" satellites selected"<<Log::endl;
    state.changedNormalEquationInfo = TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
