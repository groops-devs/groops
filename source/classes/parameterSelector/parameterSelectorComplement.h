/***********************************************/
/**
* @file parameterSelectorComplement.h
*
* @brief Parameter index vector from a complement of other parameter selector(s).
*
* @author Sebastian Strasser
* @date 2019-09-27
*
*/
/***********************************************/

#ifndef __GROOPS_PARAMETERSELECTORCOMPLEMENT__
#define __GROOPS_PARAMETERSELECTORCOMPLEMENT__

// Latex documentation
#ifdef DOCSTRING_ParameterSelector
static const char *docstringParameterSelectorComplement = R"(
\subsection{Complement}\label{parameterSelectorType:complement}
Parameter index vector from a complement of other parameter selector(s).
)";
#endif

/***********************************************/

#include "base/import.h"
#include "classes/parameterSelector/parameterSelector.h"

/***** CLASS ***********************************/

/** @brief Parameter index vector from a complement of other parameter selector(s).
* @ingroup parameterSelectorGroup
* @see ParameterSelector */
class ParameterSelectorComplement : public ParameterSelectorBase
{
  ParameterSelectorPtr parameterSelector;

public:
  ParameterSelectorComplement(Config &config);
  std::vector<UInt> indexVector(const std::vector<ParameterName> &parameterNames, VariableList varList);
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline ParameterSelectorComplement::ParameterSelectorComplement(Config &config)
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

inline std::vector<UInt> ParameterSelectorComplement::indexVector(const std::vector<ParameterName> &parameterNames, VariableList /*varList*/)
{
  try
  {
    return ParameterSelector::indexVectorComplement(parameterSelector->indexVector(parameterNames), parameterNames.size());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
