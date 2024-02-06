/***********************************************/
/**
* @file conditionOr.h
*
* @brief One of the conditions must be met (with short-circuit evaluation).
*
* @author Torsten Mayer-Guerr
* @date 2018-08-18
*
*/
/***********************************************/

#ifndef __GROOPS_CONDITIONOR__
#define __GROOPS_CONDITIONOR__

// Latex documentation
#ifdef DOCSTRING_Condition
static const char *docstringConditionOr = R"(
\subsection{Or}
One of the conditions must be met (with short-circuit evaluation).
)";
#endif

/***********************************************/

#include "base/import.h"
#include "classes/condition/condition.h"

/***** CLASS ***********************************/

/** @brief One of the conditions must be met (with short-circuit evaluation).
* @ingroup conditionGroup
* @see Condition */
class ConditionOr : public Condition
{
  std::vector<ConditionPtr> conditionPtr;

public:
  ConditionOr(Config &config);

  Bool condition(const VariableList &varList) const;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline ConditionOr::ConditionOr(Config &config)
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

inline Bool ConditionOr::condition(const VariableList &varList) const
{
  try
  {
    for(UInt i=0; i<conditionPtr.size(); i++)
      if(conditionPtr.at(i)->condition(varList))
        return TRUE;
    return FALSE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
