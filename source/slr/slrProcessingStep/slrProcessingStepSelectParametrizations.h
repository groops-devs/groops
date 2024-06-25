/***********************************************/
/**
* @file slrProcessingStepSelectParametrizations.h
*
* @brief SLR processing step: SelectParametrizations.
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#ifndef __GROOPS_SLRPROCESSINGSTEPSELECTPARAMETRIZATIONS__
#define __GROOPS_SLRPROCESSINGSTEPSELECTPARAMETRIZATIONS__

// Latex documentation
#ifdef DOCSTRING_SlrProcessingStep
static const char *docstringSlrProcessingStepSelectParametrizations = R"(
\subsection{SelectParametrizations}\label{slrProcessingStepType:selectParametrizations}
Enable/disable parameter groups and constraint groups for subsequent steps,
e.g. \configClass{processingStep:estimate}{slrProcessingStepType:estimate} or
\configClass{processingStep:writeResults}{slrProcessingStepType:writeResults}.
The \config{name} and \config{nameConstraint} of these groups
are defined in \configClass{parametrizations}{slrParametrizationType}.
Prior models or previously estimated parameters used as new apriori $\M x_0$ values are unaffected
and they are always reduced from the observations. This means all unselected parameters are kept fixed
to their last result.
)";
#endif

/***********************************************/

#include "config/config.h"
#include "base/string.h"
#include "slr/slrProcessingStep/slrProcessingStep.h"

/***** CLASS ***********************************/

/** @brief SLR processing step: SelectParametrizations.
* @ingroup slrProcessingStepGroup
* @see SlrProcessingStep */
class SlrProcessingStepSelectParametrizations : public SlrProcessingStepBase
{
public:
  class EnableDisable
  {
  public:
    Bool enable;
    std::vector<std::string> wildcards;
  };

  std::vector<EnableDisable> patterns;

  SlrProcessingStepSelectParametrizations(Config &config);
  void process(SlrProcessingStep::State &state) override;
  Bool expectInitializedParameters() const override {return FALSE;}
};

/***********************************************/

template<> inline Bool readConfig(Config &config, const std::string &name, SlrProcessingStepSelectParametrizations::EnableDisable &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  try
  {
    std::string choice;
    if(!readConfigChoice(config, name, choice, mustSet, defaultValue, annotation))
      return FALSE;
    if(readConfigChoiceElement(config, "enable", choice, ""))
    {
      var.enable = TRUE;
      readConfig(config, "name", var.wildcards, Config::MUSTSET, "*", "wildcards: * and ?");
    }
    if(readConfigChoiceElement(config, "disable", choice, ""))
    {
      var.enable = FALSE;
      readConfig(config, "name", var.wildcards, Config::MUSTSET, "*", "wildcards: * and ?");
    }
    endChoice(config);
    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline SlrProcessingStepSelectParametrizations::SlrProcessingStepSelectParametrizations(Config &config)
{
  try
  {
    readConfig(config, "parametrization", patterns, Config::MUSTSET, "", "");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void SlrProcessingStepSelectParametrizations::process(SlrProcessingStep::State &state)
{
  try
  {
    logStatus<<"=== select parametrization ======================================="<<Log::endl;
    state.changedNormalEquationInfo = TRUE;
    for(const auto &pattern : patterns)
      for(const auto &wildcard: pattern.wildcards)
      {
        state.normalEquationInfo.enableParametrizations.emplace_back(String::wildcard2regex(wildcard), pattern.enable);
        logInfo<<"  "<<(pattern.enable ? "enable " : "disable")<<" parametrization: "<<wildcard<<Log::endl;
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
