/***********************************************/
/**
* @file conditionExpression.h
*
* @brief Evaluate expression.
*
* @author Torsten Mayer-Guerr
* @date 2018-08-18
*
*/
/***********************************************/

#ifndef __GROOPS_CONDITIONEXPRESSION__
#define __GROOPS_CONDITIONEXPRESSION__

// Latex documentation
#ifdef DOCSTRING_Condition
static const char *docstringConditionExpression = R"(
\subsection{Expression}
Evaluate expression.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "classes/condition/condition.h"

/***** CLASS ***********************************/

/** @brief Evaluate expression.
* @ingroup conditionGroup
* @see Condition */
class ConditionExpression : public Condition
{
  ExpressionVariablePtr expr;

public:
  ConditionExpression(Config &config);

  Bool condition(const VariableList &varList) const;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline ConditionExpression::ConditionExpression(Config &config)
{
  try
  {
    readConfig(config, "expression", expr, Config::MUSTSET, "",  "");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Bool ConditionExpression::condition(const VariableList &varList) const
{
  try
  {
    return (expr->evaluate(varList) != 0.);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
