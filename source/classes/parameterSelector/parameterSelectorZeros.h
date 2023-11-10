/***********************************************/
/**
* @file parameterSelectorZeros.h
*
* @brief Expand parameter index vector by adding zero parameters.
*
* @author Sebastian Strasser
* @date 2018-05-08
*
*/
/***********************************************/

#ifndef __GROOPS_PARAMETERSELECTORZEROS__
#define __GROOPS_PARAMETERSELECTORZEROS__

// Latex documentation
#ifdef DOCSTRING_ParameterSelector
static const char *docstringParameterSelectorZeros = R"(
\subsection{Zeros}
Expand parameter index vector by adding zero parameters.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "classes/parameterSelector/parameterSelector.h"

/***** CLASS ***********************************/

/** @brief Expand parameter index vector by adding zero parameters.
* @ingroup parameterSelectorGroup
* @see ParameterSelector */
class ParameterSelectorZeros : public ParameterSelectorBase
{
  ExpressionVariablePtr exprCount;

public:
  ParameterSelectorZeros(Config &config);
  std::vector<UInt> indexVector(const std::vector<ParameterName> &parameterNames);
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline ParameterSelectorZeros::ParameterSelectorZeros(Config &config)
{
  try
  {
    readConfig(config, "count", exprCount, Config::MUSTSET, "", "count of zero parameters (variables: length)");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline std::vector<UInt> ParameterSelectorZeros::indexVector(const std::vector<ParameterName> &parameterNames)
{
  try
  {
    VariableList varList;
    varList.setVariable("length", static_cast<Double>(parameterNames.size()));
    const UInt count = static_cast<UInt>(exprCount->evaluate(varList));
    return std::vector<UInt>(count, NULLINDEX);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
