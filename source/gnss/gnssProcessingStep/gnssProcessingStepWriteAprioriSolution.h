/***********************************************/
/**
* @file gnssProcessingStepWriteAprioriSolution.h
*
* @brief GNSS processing step: WriteAprioriSolution.
*
* @author Torsten Mayer-Guerr
* @date 2021-09-05
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPROCESSINGSTEPWRITEAPRIORISOLUTION__
#define __GROOPS_GNSSPROCESSINGSTEPWRITEAPRIORISOLUTION__

// Latex documentation
#ifdef DOCSTRING_GnssProcessingStep
static const char *docstringGnssProcessingStepWriteAprioriSolution = R"(
\subsection{WriteAprioriSolution}\label{gnssProcessingStepType:writeAprioriSolution}
Writes the current apriori vector $\M x_0$
(see \configClass{parametrizations}{gnssParametrizationType}).
If \configClass{remainingParameters}{parameterSelectorType}
is set only the selected parameters are written.
)";
#endif

/***********************************************/

#include "config/config.h"
#include "files/fileMatrix.h"
#include "files/fileParameterName.h"
#include "classes/parameterSelector/parameterSelector.h"
#include "gnss/gnssProcessingStep/gnssProcessingStep.h"

/***** CLASS ***********************************/

/** @brief GNSS processing step: WriteAprioriSolution.
* @ingroup gnssProcessingStepGroup
* @see GnssProcessingStep */
class GnssProcessingStepWriteAprioriSolution : public GnssProcessingStepBase
{
  FileName             fileNameApriori, fileNameParameterNames;
  ParameterSelectorPtr parameterSelector;

public:
  GnssProcessingStepWriteAprioriSolution(Config &config);
  void process(GnssProcessingStep::State &state) override;
};

/***********************************************/

inline GnssProcessingStepWriteAprioriSolution::GnssProcessingStepWriteAprioriSolution(Config &config)
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

inline void GnssProcessingStepWriteAprioriSolution::process(GnssProcessingStep::State &state)
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
      const Vector x0 = reorder(state.gnss->aprioriParameter(state.normalEquationInfo), indexVector);
      if(Parallel::isMaster(state.normalEquationInfo.comm))
        writeFileMatrix(fileNameApriori, x0);
    }

    if(!fileNameParameterNames.empty() && Parallel::isMaster(state.normalEquationInfo.comm))
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
