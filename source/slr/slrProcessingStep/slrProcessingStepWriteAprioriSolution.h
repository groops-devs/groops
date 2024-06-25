/***********************************************/
/**
* @file slrProcessingStepWriteAprioriSolution.h
*
* @brief SLR processing step: WriteAprioriSolution.
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#ifndef __GROOPS_SLRPROCESSINGSTEPWRITEAPRIORISOLUTION__
#define __GROOPS_SLRPROCESSINGSTEPWRITEAPRIORISOLUTION__

// Latex documentation
#ifdef DOCSTRING_SlrProcessingStep
static const char *docstringSlrProcessingStepWriteAprioriSolution = R"(
\subsection{WriteAprioriSolution}\label{slrProcessingStepType:writeAprioriSolution}
Writes the current apriori vector $\M x_0$
(see \configClass{parametrizations}{slrParametrizationType}).
If \configClass{remainingParameters}{parameterSelectorType}
is set only the selected parameters are written.
)";
#endif

/***********************************************/

#include "config/config.h"
#include "files/fileMatrix.h"
#include "files/fileParameterName.h"
#include "classes/parameterSelector/parameterSelector.h"
#include "slr/slrProcessingStep/slrProcessingStep.h"

/***** CLASS ***********************************/

/** @brief SLR processing step: WriteAprioriSolution.
* @ingroup slrProcessingStepGroup
* @see SlrProcessingStep */
class SlrProcessingStepWriteAprioriSolution : public SlrProcessingStepBase
{
  FileName             fileNameApriori, fileNameParameterNames;
  ParameterSelectorPtr parameterSelector;

public:
  SlrProcessingStepWriteAprioriSolution(Config &config);
  void process(SlrProcessingStep::State &state) override;
};

/***********************************************/

inline SlrProcessingStepWriteAprioriSolution::SlrProcessingStepWriteAprioriSolution(Config &config)
{
  try
  {
    readConfig(config, "outputfileAprioriSolution", fileNameApriori,        Config::OPTIONAL, "output/x0_{loopTime:%D}.txt",            "a priori parameters");
    readConfig(config, "outputfileParameterNames",  fileNameParameterNames, Config::OPTIONAL, "output/parameterNames_{loopTime:%D}.txt", "parameter names");
    readConfig(config, "remainingParameters",       parameterSelector,      Config::OPTIONAL, "",  "parameter order/selection of output normal equations");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void SlrProcessingStepWriteAprioriSolution::process(SlrProcessingStep::State &state)
{
  try
  {
    std::vector<UInt> indexVector(state.normalEquationInfo.parameterCount());
    if(parameterSelector)
      indexVector = parameterSelector->indexVector(state.normalEquationInfo.parameterNames());
    else
      std::iota(indexVector.begin(), indexVector.end(), 0);

    if(!fileNameApriori.empty())
    {
      logStatus<<"write apriori solution to <"<<fileNameApriori<<">"<<Log::endl;
      const Vector x0 = reorder(state.slr->aprioriParameter(state.normalEquationInfo), indexVector);
      writeFileMatrix(fileNameApriori, x0);
    }

    if(!fileNameParameterNames.empty())
    {
      logStatus<<"write parameter name file <"<<fileNameParameterNames<<">"<<Log::endl;
      std::vector<ParameterName> parameterNames;
      for(UInt i=0; i<indexVector.size(); i++)
        parameterNames.push_back((indexVector.at(i) != NULLINDEX) ? state.normalEquationInfo.parameterNames().at(indexVector.at(i)) : ParameterName());
      writeFileParameterName(fileNameParameterNames, parameterNames);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
