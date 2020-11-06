/***********************************************/
/**
* @file parameterNamesSelection.h
*
* @brief Selection of parameter names.
*
* @author Torsten Mayer-Guerr
* @date 2020-05-29
*/
/***********************************************/

#ifndef __GROOPS_PARAMETERNAMESSELECTION__
#define __GROOPS_PARAMETERNAMESSELECTION__

// Latex documentation
static const char *docstringParameterNamesSelection = R"(
\subsection{Selection}
Select a subset of \configClass{parameterName}{parameterNamesType}s
using \configClass{parameterSelection}{parameterSelectorType}.
)";

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "classes/parameterSelector/parameterSelector.h"
#include "classes/parameterNames/parameterNames.h"

/***** CLASS ***********************************/

/** @brief Selection of parameter names.
* @ingroup parameterNamesGroup
* @see ParameterNames */
class ParameterNamesSelection : public ParameterNamesBase
{
public:
  ParameterNamesSelection(Config &config)
  {
    try
    {
      ParameterNamesPtr    parameterNames;
      ParameterSelectorPtr selector;

      readConfig(config, "parameterName",      parameterNames, Config::MUSTSET, "", "");
      readConfig(config, "parameterSelection", selector,       Config::MUSTSET, "", "parameter order/selection");
      if(isCreateSchema(config)) return;

      for(UInt i : selector->indexVector(parameterNames->parameterNames()))
        names.push_back((i != NULLINDEX) ? parameterNames->parameterNames().at(i) : ParameterName());
    }
    catch(std::exception &e)
    {
      GROOPS_RETHROW(e)
    }
  }
};

/***********************************************/

#endif /* __GROOPS__ */
