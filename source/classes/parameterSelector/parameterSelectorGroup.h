/***********************************************/
/**
* @file parameterSelectorGroup.h
*
* @brief Groups a set parameter selectors.
*
* @author Torsten Mayer-Guerr
* @date 2026-93-22
*
*/
/***********************************************/

#ifndef __GROOPS_PARAMETERSELECTORGROUP__
#define __GROOPS_PARAMETERSELECTORGROUP__

// Latex documentation
#ifdef DOCSTRING_ParameterSelector
static const char *docstringParameterSelectorGroup = R"(
\subsection{Group}\label{parameterSelectorType:group}
Groups a set of \configClass{parameterSelection}{parameterSelectorType}s and has no further effect itself.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "classes/parameterSelector/parameterSelector.h"

/***** CLASS ***********************************/

/** @brief Groups a set parameter selectors.
* @ingroup parameterSelectorGroup
* @see ParameterSelector */
class ParameterSelectorGroup : public ParameterSelectorBase
{
  ParameterSelectorPtr parameterSelector;

public:
  ParameterSelectorGroup(Config &config);
  std::vector<UInt> indexVector(const std::vector<ParameterName> &parameterNames);
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline ParameterSelectorGroup::ParameterSelectorGroup(Config &config)
{
  try
  {
    readConfig(config, "parameterSelection", parameterSelector, Config::MUSTSET, "", "parameter order/selection");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline std::vector<UInt> ParameterSelectorGroup::indexVector(const std::vector<ParameterName> &parameterNames)
{
  try
  {
    return parameterSelector->indexVector(parameterNames);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
