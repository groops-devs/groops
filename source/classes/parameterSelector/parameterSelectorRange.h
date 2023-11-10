/***********************************************/
/**
* @file parameterSelectorRange.h
*
* @brief Parameter index vector from range.
*
* @author Sebastian Strasser
* @date 2018-05-08
*
*/
/***********************************************/

#ifndef __GROOPS_PARAMETERSELECTORRANGE__
#define __GROOPS_PARAMETERSELECTORRANGE__

// Latex documentation
#ifdef DOCSTRING_ParameterSelector
static const char *docstringParameterSelectorRange = R"(
\subsection{Range}
Parameter index vector from range.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "classes/parameterSelector/parameterSelector.h"

/***** CLASS ***********************************/

/** @brief Parameter index vector from range.
* @ingroup parameterSelectorGroup
* @see ParameterSelector */
class ParameterSelectorRange : public ParameterSelectorBase
{
  ExpressionVariablePtr exprStart, exprCount;

public:
  ParameterSelectorRange(Config &config);
  std::vector<UInt> indexVector(const std::vector<ParameterName> &parameterNames);
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline ParameterSelectorRange::ParameterSelectorRange(Config &config)
{
  try
  {
    readConfig(config, "start",   exprStart, Config::MUSTSET,  "0", "start at this index (variables: length)");
    readConfig(config, "count",   exprCount, Config::OPTIONAL, "",  "count of parameters, default: all parameters to the end (variables: length)");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline std::vector<UInt> ParameterSelectorRange::indexVector(const std::vector<ParameterName> &parameterNames)
{
  try
  {
    const UInt parameterCount = parameterNames.size();
    VariableList varList;
    varList.setVariable("length", static_cast<Double>(parameterCount));
    const UInt start = static_cast<UInt>(exprStart->evaluate(varList));
    const UInt count = (exprCount ? static_cast<UInt>(exprCount->evaluate(varList)) : parameterCount-start);

    if(start >= parameterCount)
      return {};
    if(start+count > parameterCount)
      throw(Exception("range exceeds parameter count: "+(start+count)%"%i >= "s+parameterCount%"%i"s));

    std::vector<UInt> vector(count);
    std::iota(vector.begin(), vector.end(), start);

    return vector;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
