/***********************************************/
/**
* @file slrProcessingStepWriteResults.h
*
* @brief SLR processing step: WriteResults.
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#ifndef __GROOPS_SLRPROCESSINGSTEPWRITERESULTS__
#define __GROOPS_SLRPROCESSINGSTEPWRITERESULTS__

// Latex documentation
#ifdef DOCSTRING_SlrProcessingStep
static const char *docstringSlrProcessingStepWriteResults = R"(
\subsection{WriteResults}\label{slrProcessingStepType:writeResults}
In this step all \config{outputfiles} defined in \configClass{parametrizations}{slrParametrizationType}
are written. It considers the settings of
\configClass{processingStep:selectParametrizations}{slrProcessingStepType:selectParametrizations}
and \configClass{processingStep:selectStations}{slrProcessingStepType:selectStations}.

It is usually the last processing step, but can also be used at other points in the
processing in combination with \config{suffix} to write intermediate results.
)";
#endif

/***********************************************/

#include "config/config.h"
#include "slr/slrProcessingStep/slrProcessingStep.h"

/***** CLASS ***********************************/

/** @brief SLR processing step: WriteResults.
* @ingroup slrProcessingStepGroup
* @see SlrProcessingStep */
class SlrProcessingStepWriteResults : public SlrProcessingStepBase
{
  std::string suffix;

public:
  SlrProcessingStepWriteResults(Config &config);
  void process(SlrProcessingStep::State &state) override;
};

/***********************************************/

inline SlrProcessingStepWriteResults::SlrProcessingStepWriteResults(Config &config)
{
  readConfig(config, "suffix", suffix, Config::OPTIONAL, "", "appended to every output file name (e.g. orbit.G01.suffix.dat)");
}

/***********************************************/

inline void SlrProcessingStepWriteResults::process(SlrProcessingStep::State &state)
{
  try
  {
    logStatus<<"=== write results  =========================================="<<Log::endl;
    state.slr->writeResults(state.normalEquationInfo, !suffix.empty() ? "."+suffix : "");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
