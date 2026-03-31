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
Supports wildcards * for any number of characters and ? for exactly one character.
Files smaller than \config{minSize} are treated as non-existent.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "parallel/parallel.h"
#include "inputOutput/system.h"
#include "classes/condition/condition.h"

/***** CLASS ***********************************/

/** @brief Check for a file or directory existing.
* @ingroup conditionGroup
* @see Condition */
class ConditionFileExist : public Condition
{
  FileName fileName;
  UInt     minSize;

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
    minSize = 0;

    readConfig(config, "file",        fileName, Config::MUSTSET,  "", "supports wildcards: * and ?");
    readConfig(config, "minimumSize", minSize,  Config::OPTIONAL, "", "minimum file size in byte.");
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
    const std::vector<FileName> fileList = System::fileList(fileName(varList));
    for(auto entry : fileList)
      if(System::isDirectory(entry) || (System::fileSize(entry) >= minSize))
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
