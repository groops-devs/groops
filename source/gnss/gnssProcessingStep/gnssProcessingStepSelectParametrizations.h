/***********************************************/
/**
* @file gnssProcessingStepSelectParametrizations.h
*
* @brief GNSS processing step: SelectParametrizations.
*
* @author Torsten Mayer-Guerr
* @date 2021-09-05
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPROCESSINGSTEPSELECTPARAMETRIZATIONS__
#define __GROOPS_GNSSPROCESSINGSTEPSELECTPARAMETRIZATIONS__

// Latex documentation
#ifdef DOCSTRING_GnssProcessingStep
static const char *docstringGnssProcessingStepSelectParametrizations = R"(
\subsection{SelectParametrizations}\label{gnssProcessingStepType:selectParametrizations}
Enable/disable parameter groups and constraint groups for subsequent steps,
e.g. \configClass{processingStep:estimate}{gnssProcessingStepType:estimate} or
\configClass{processingStep:writeResults}{gnssProcessingStepType:writeResults}.
The \config{name} and \config{nameConstraint} of these groups
are defined in \configClass{parametrizations}{gnssParametrizationType}.
Prior models or previously estimated parameters used as new apriori $\M x_0$ values are unaffected
and they are always reduced from the observations. This means all unselected parameters are kept fixed
to their last result.

An example would be to process at a 5-minute sampling using
\configClass{processingStep:selectEpochs}{gnssProcessingStepType:selectEpochs}
and then at the end to densify the clock parameters to the full 30-second observation sampling
while keeping all other parameters fixed
(\config{disable}=\verb|*|, \config{enable}=\verb|*.clock*|, \config{enable}=\verb|parameter.STEC|).
)";
#endif

/***********************************************/

#include "config/config.h"
#include "base/string.h"
#include "gnss/gnssProcessingStep/gnssProcessingStep.h"

/***** CLASS ***********************************/

/** @brief GNSS processing step: SelectParametrizations.
* @ingroup gnssProcessingStepGroup
* @see GnssProcessingStep */
class GnssProcessingStepSelectParametrizations : public GnssProcessingStepBase
{
public:
  class EnableDisable
  {
  public:
    Bool enable;
    std::vector<std::string> wildcards;
  };

  std::vector<EnableDisable> patterns;

  GnssProcessingStepSelectParametrizations(Config &config);
  void process(GnssProcessingStep::State &state) override;
  Bool expectInitializedParameters() const override {return FALSE;}
};

/***********************************************/

template<> inline Bool readConfig(Config &config, const std::string &name, GnssProcessingStepSelectParametrizations::EnableDisable &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
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

inline GnssProcessingStepSelectParametrizations::GnssProcessingStepSelectParametrizations(Config &config)
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

void GnssProcessingStepSelectParametrizations::process(GnssProcessingStep::State &state)
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
