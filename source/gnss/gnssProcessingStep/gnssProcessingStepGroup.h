/***********************************************/
/**
* @file gnssProcessingStepGroup.h
*
* @brief GNSS processing step: Group.
*
* @author Torsten Mayer-Guerr
* @date 2021-09-05
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPROCESSINGSTEPGROUP__
#define __GROOPS_GNSSPROCESSINGSTEPGROUP__

// Latex documentation
#ifdef DOCSTRING_GnssProcessingStep
static const char *docstringGnssProcessingStepGroup = R"(
\subsection{Group}\label{gnssProcessingStepType:group}
Perform these processing steps. This step can be used to structure complex processing flows.
The \configClass{select..}{gnssProcessingStepType:selectParametrizations} processing steps
defined within a group only affect the steps within this group.
)";
#endif

/***********************************************/

#include "config/config.h"
#include "gnss/gnssProcessingStep/gnssProcessingStep.h"

/***** CLASS ***********************************/

/** @brief GNSS processing step: Group.
* @ingroup gnssProcessingStepGroup
* @see GnssProcessingStep */
class GnssProcessingStepGroup : public GnssProcessingStepBase
{
  GnssProcessingStepPtr processingSteps;

public:
  GnssProcessingStepGroup(Config &config);
  void process(GnssProcessingStep::State &state) override;
  Bool expectInitializedParameters() const override {return FALSE;}
};

/***********************************************/

inline GnssProcessingStepGroup::GnssProcessingStepGroup(Config &config)
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

void GnssProcessingStepGroup::process(GnssProcessingStep::State &state)
{
  try
  {
    logStatus<<"=== process group ==========================================="<<Log::endl;
    GnssNormalEquationInfo normalEquationInfoOld = state.normalEquationInfo;
    processingSteps->process(state);
    std::swap(state.normalEquationInfo, normalEquationInfoOld);
    state.changedNormalEquationInfo = TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
