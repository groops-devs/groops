/***********************************************/
/**
* @file parameterSelectorNames.h
*
* @brief Parameter index vector from list of parameter names.
*
* @author Sebastian Strasser
* @date 2018-05-08
*
*/
/***********************************************/

#ifndef __GROOPS_PARAMETERSELECTORMANUALLIST__
#define __GROOPS_PARAMETERSELECTORMANUALLIST__

// Latex documentation
#ifdef DOCSTRING_ParameterSelector
static const char *docstringParameterSelectorNames = R"(
\subsection{Names}
Parameter index vector from list of parameter names.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "classes/parameterNames/parameterNames.h"
#include "classes/parameterSelector/parameterSelector.h"

/***** CLASS ***********************************/

/** @brief Parameter index vector from list of parameter names.
* @ingroup parameterSelectorGroup
* @see ParameterSelector */
class ParameterSelectorNames : public ParameterSelectorBase
{
  std::vector<ParameterName> requestedNames;

public:
  ParameterSelectorNames(Config &config);
  std::vector<UInt> indexVector(const std::vector<ParameterName> &ParameterNames, VariableList varList);
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline ParameterSelectorNames::ParameterSelectorNames(Config &config)
{
  try
  {
    ParameterNamesPtr parameterNames;
    readConfig(config, "parameterName", parameterNames, Config::MUSTSET,  "",  "");
    if(isCreateSchema(config)) return;

    requestedNames = parameterNames->parameterNames();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline std::vector<UInt> ParameterSelectorNames::indexVector(const std::vector<ParameterName> &parameterNames, VariableList /*varList*/)
{
  try
  {
    std::vector<UInt> vector;
    for(const auto &name : requestedNames)
    {
      auto iter = std::find(parameterNames.begin(), parameterNames.end(), name);
      if(iter == parameterNames.end())
        vector.push_back(NULLINDEX);
      else
        vector.push_back(std::distance(parameterNames.begin(), iter));
    }

    return vector;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
