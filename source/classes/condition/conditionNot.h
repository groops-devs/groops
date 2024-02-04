/***********************************************/
/**
* @file conditionNot.h
*
* @brief The result of the condition is inverted.
*
* @author Torsten Mayer-Guerr
* @date 2019-03-20
*
*/
/***********************************************/

#ifndef __GROOPS_CONDITIONNOT__
#define __GROOPS_CONDITIONNOT__

// Latex documentation
#ifdef DOCSTRING_Condition
static const char *docstringConditionNot = R"(
\subsection{Not}
The result of the condition is inverted.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "classes/condition/condition.h"

/***** CLASS ***********************************/

/** @brief The result of the condition is inverted.
* @ingroup conditionGroup
* @see Condition */
class ConditionNot : public Condition
{
  ConditionPtr conditionPtr;

public:
  ConditionNot(Config &config);

  Bool condition(const VariableList &varList) const;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline ConditionNot::ConditionNot(Config &config)
{
  try
  {
    readConfig(config, "condition", conditionPtr, Config::MUSTSET, "", "");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Bool ConditionNot::condition(const VariableList &varList) const
{
  try
  {
    return !conditionPtr->condition(varList);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
