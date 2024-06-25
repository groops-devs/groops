/***********************************************/
/**
* @file slrProcessingStepGroup.h
*
* @brief SLR processing step: Group.
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#ifndef __GROOPS_SLRPROCESSINGSTEPGROUP__
#define __GROOPS_SLRPROCESSINGSTEPGROUP__

// Latex documentation
#ifdef DOCSTRING_SlrProcessingStep
static const char *docstringSlrProcessingStepGroup = R"(
\subsection{Group}\label{slrProcessingStepType:group}
Perform these processing steps. This step can be used to structure complex processing flows.
The \configClass{select..}{slrProcessingStepType:selectParametrizations} processing steps
defined within a group only affect the steps within this group.
)";
#endif

/***********************************************/

#include "config/config.h"
#include "slr/slrProcessingStep/slrProcessingStep.h"

/***** CLASS ***********************************/

/** @brief SLR processing step: Group.
* @ingroup slrProcessingStepGroup
* @see SlrProcessingStep */
class SlrProcessingStepGroup : public SlrProcessingStepBase
{
  SlrProcessingStepPtr processingSteps;

public:
  SlrProcessingStepGroup(Config &config);
  void process(SlrProcessingStep::State &state) override;
  Bool expectInitializedParameters() const override {return FALSE;}
};

/***********************************************/

inline SlrProcessingStepGroup::SlrProcessingStepGroup(Config &config)
{
  try
  {
    readConfig(config, "processingStep", processingSteps, Config::MUSTSET, "", "steps are processed consecutively");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void SlrProcessingStepGroup::process(SlrProcessingStep::State &state)
{
  try
  {
    logStatus<<"=== process group ==========================================="<<Log::endl;
    SlrNormalEquationInfo normalEquationInfoOld = state.normalEquationInfo;
    processingSteps->process(state);
    state.normalEquationInfo = std::move(normalEquationInfoOld);
    state.changedNormalEquationInfo = TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
