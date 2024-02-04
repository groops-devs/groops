/***********************************************/
/**
* @file conditionAnd.h
*
* @brief All conditions must be met (with short-circuit evaluation).
*
* @author Torsten Mayer-Guerr
* @date 2018-08-18
*
*/
/***********************************************/

#ifndef __GROOPS_CONDITIONAND__
#define __GROOPS_CONDITIONAND__

// Latex documentation
#ifdef DOCSTRING_Condition
static const char *docstringConditionAnd = R"(
\subsection{And}
All conditions must be met (with short-circuit evaluation).
)";
#endif

/***********************************************/

#include "base/import.h"
#include "classes/condition/condition.h"

/***** CLASS ***********************************/

/** @brief All conditions must be met (with short-circuit evaluation).
* @ingroup conditionGroup
* @see Condition */
class ConditionAnd : public Condition
{
  std::vector<ConditionPtr> conditionPtr;

public:
  ConditionAnd(Config &config);

  Bool condition(const VariableList &varList) const;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline ConditionAnd::ConditionAnd(Config &config)
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

inline Bool ConditionAnd::condition(const VariableList &varList) const
{
  try
  {
    for(UInt i=0; i<conditionPtr.size(); i++)
      if(!conditionPtr.at(i)->condition(varList))
        return FALSE;
    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
