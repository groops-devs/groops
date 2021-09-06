/***********************************************/
/**
* @file gnssProcessingStepWriteResults.h
*
* @brief GNSS processing step: WriteResults.
*
* @author Torsten Mayer-Guerr
* @date 2021-09-05
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPROCESSINGSTEPWRITERESULTS__
#define __GROOPS_GNSSPROCESSINGSTEPWRITERESULTS__

// Latex documentation
#ifdef DOCSTRING_GnssProcessingStep
static const char *docstringGnssProcessingStepWriteResults = R"(
\subsection{WriteResults}\label{gnssProcessingStepType:writeResults}
In this step all \config{outputfiles} defined in \configClass{parametrizations}{gnssParametrizationType}
are written. It considers the settings of
\configClass{processingStep:selectParametrizations}{gnssProcessingStepType:selectParametrizations},
\configClass{processingStep:selectEpochs}{gnssProcessingStepType:selectEpochs}, and
\configClass{processingStep:selectReceivers}{gnssProcessingStepType:selectReceivers}.

It is usually the last processing step, but can also be used at other points in the
processing in combination with \config{suffix} to write intermediate results, for example
before \configClass{gnssProcessingStep:resolveAmbiguities}{gnssProcessingStepType:resolveAmbiguities} to
output the float solution.
)";
#endif

/***********************************************/

#include "config/config.h"
#include "gnss/gnssProcessingStep/gnssProcessingStep.h"

/***** CLASS ***********************************/

/** @brief GNSS processing step: WriteResults.
* @ingroup gnssProcessingStepGroup
* @see GnssProcessingStep */
class GnssProcessingStepWriteResults : public GnssProcessingStepBase
{
  std::string suffix;

public:
  GnssProcessingStepWriteResults(Config &config);
  void process(GnssProcessingStep::State &state) override;
};

/***********************************************/

inline GnssProcessingStepWriteResults::GnssProcessingStepWriteResults(Config &config)
{
  readConfig(config, "suffix", suffix, Config::OPTIONAL, "", "appended to every output file name (e.g. orbit.G01.suffix.dat)");
}

/***********************************************/

void GnssProcessingStepWriteResults::process(GnssProcessingStep::State &state)
{
  try
  {
    logStatus<<"=== write results  =========================================="<<Log::endl;
    if(state.changedNormalEquationInfo)
      state.gnss->initParameter(state.normalEquationInfo);
    state.changedNormalEquationInfo = FALSE;
    state.gnss->writeResults(state.normalEquationInfo, !suffix.empty() ? "."+suffix : "");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
