/***********************************************/
/**
* @file conditionFileExist.h
*
* @brief Check for a file or directory existing.
*
* @author Torsten Mayer-Guerr
* @date 2018-08-18
*
*/
/***********************************************/

#ifndef __GROOPS_CONDITIONFILEEXIST__
#define __GROOPS_CONDITIONFILEEXIST__

// Latex documentation
#ifdef DOCSTRING_Condition
static const char *docstringConditionFileExist = R"(
\subsection{FileExist}
Check for a file or directory existing.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "parallel/parallel.h"
#include "parser/stringParser.h"
#include "inputOutput/system.h"
#include "classes/condition/condition.h"

/***** CLASS ***********************************/

/** @brief Check for a file or directory existing.
* @ingroup ConditionGroup
* @see Condition */
class ConditionFileExist : public Condition
{
  FileName fileName;

public:
  ConditionFileExist(Config &config);

  Bool condition(const VariableList &varList) const;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline ConditionFileExist::ConditionFileExist(Config &config)
{
  try
  {
    readConfig(config, "file", fileName, Config::MUSTSET,  "",  "");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Bool ConditionFileExist::condition(const VariableList &varList) const
{
  try
  {
    Bool result;
    if(Parallel::isMaster())
      result = System::exists(StringParser::parse(fileName, varList));
    Parallel::broadCast(result);
    return result;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
