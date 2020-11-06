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
  std::string nameString, nameIndex, nameCount;
  std::vector<std::string> strings;

public:
  LoopManualList(Config &config);

  UInt count() const override;
  void init(VariableList &varList) override;
  void setValues(VariableList &varList) override;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline LoopManualList::LoopManualList(Config &config) : Loop()
{
  try
  {
    readConfig(config, "string",             strings,    Config::MUSTSET,   "",           "explicit list of strings");
    readConfig(config, "variableLoopString", nameString, Config::OPTIONAL,  "loopString", "name of the variable to be replaced");
    readConfig(config, "variableLoopIndex",  nameIndex,  Config::OPTIONAL,  "",           "variable with index of current iteration (starts with zero)");
    readConfig(config, "variableLoopCount",  nameCount,  Config::OPTIONAL,  "",           "variable with total number of iterations");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline UInt LoopManualList::count() const
{
  return strings.size();
}

/***********************************************/

inline void LoopManualList::init(VariableList &varList)
{
  if(index == NULLINDEX)
  {
    if(!nameString.empty()) addVariable(nameString, varList);
    if(!nameIndex.empty())  addVariable(nameIndex,  varList);
    if(!nameCount.empty())  addVariable(nameCount,  varList);
  }
  index = 0;
}

/***********************************************/

inline void LoopManualList::setValues(VariableList &varList)
{
  if(!nameString.empty()) varList[nameString]->setValue(strings.at(index));
  if(!nameIndex.empty())  varList[nameIndex]->setValue(index);
  if(!nameCount.empty())  varList[nameCount]->setValue(count());
}

/***********************************************/


#endif
