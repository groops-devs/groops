/***********************************************/
/**
* @file loopManualList.h
*
* @brief Loop over list of strings.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2017-01-27
*
*/
/***********************************************/

#ifndef __GROOPS_LOOPMANUALLIST__
#define __GROOPS_LOOPMANUALLIST__

// Latex documentation
#ifdef DOCSTRING_Loop
static const char *docstringLoopManualList = R"(
\subsection{ManualList}\label{loopType:manualList}
Loop over list of strings.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "classes/loop/loop.h"

/***** CLASS ***********************************/

/** @brief Loop over list of strings.
* @ingroup LoopGroup
* @see Loop */
class LoopManualList : public Loop
{
  std::string              nameString, nameIndex, nameCount;
  std::vector<std::string> strings;
  UInt                     index;

public:
  LoopManualList(Config &config);

  UInt count() const override {return strings.size();}
  Bool iteration(VariableList &varList) override;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline LoopManualList::LoopManualList(Config &config)
{
  try
  {
    readConfig(config, "string",             strings,    Config::MUSTSET,   "",           "explicit list of strings");
    readConfig(config, "variableLoopString", nameString, Config::OPTIONAL,  "loopString", "name of the variable to be replaced");
    readConfig(config, "variableLoopIndex",  nameIndex,  Config::OPTIONAL,  "",           "variable with index of current iteration (starts with zero)");
    readConfig(config, "variableLoopCount",  nameCount,  Config::OPTIONAL,  "",           "variable with total number of iterations");
    if(isCreateSchema(config)) return;

    index = 0;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Bool LoopManualList::iteration(VariableList &varList)
{
  if(index >= count())
    return FALSE;

  if(!nameString.empty()) addVariable(nameString, strings.at(index), varList);
  if(!nameIndex.empty())  addVariable(nameIndex,  index,             varList);
  if(!nameCount.empty())  addVariable(nameCount,  count(),           varList);

  index++;
  return TRUE;
}

/***********************************************/

#endif
