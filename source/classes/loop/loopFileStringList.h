/***********************************************/
/**
* @file loopFileStringList.h
*
* @brief Loop over list of strings from files.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2017-01-27
*
*/
/***********************************************/

#ifndef __GROOPS_LOOPFILESTRINGLIST__
#define __GROOPS_LOOPFILESTRINGLIST__

// Latex documentation
#ifdef DOCSTRING_Loop
static const char *docstringLoopFileStringList = R"(
\subsection{FileStringList}\label{loopType:fileStringList}
Loop over list of strings from \file{files}{stringList}.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "files/fileStringTable.h"
#include "classes/loop/loop.h"
#include <unordered_set>

/***** CLASS ***********************************/

/** @brief Loop over list of strings from files.
* @ingroup loopGroup
* @see Loop */
class LoopFileStringList : public Loop
{
  std::string nameString, nameIndex, nameCount;
  std::vector<std::string> strings;

public:
  LoopFileStringList(Config &config);

  UInt count() const override {return strings.size();}
  Bool iteration(VariableList &varList) override;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline LoopFileStringList::LoopFileStringList(Config &config)
{
  try
  {
    std::vector<FileName> fileName;

    readConfig(config, "inputfile",          fileName,   Config::MUSTSET,  "",           "string list file");
    readConfig(config, "variableLoopString", nameString, Config::OPTIONAL, "loopString", "name of the variable to be replaced");
    readConfig(config, "variableLoopIndex",  nameIndex,  Config::OPTIONAL, "",           "variable with index of current iteration (starts with zero)");
    readConfig(config, "variableLoopCount",  nameCount,  Config::OPTIONAL, "",           "variable with total number of iterations");
    readConfigCondition(config);
    if(isCreateSchema(config)) return;

    for(UInt i=0; i<fileName.size(); i++)
    {
      std::vector<std::string> str;
      readFileStringList(fileName.at(i), str);
      strings.insert(strings.end(), str.begin(), str.end());
    }

    // DEPRECTATED
    // -----------
    // read deprecated config elements (as after if(isCreateSchema(config)) not shown in the GUI and docu)
    UInt startIndex;
    UInt countElements = MAX_UINT;
    Bool sort;
    Bool removeDuplicates;
    readConfig(config, "sort",               sort,             Config::DEFAULT,  "0", "sort entries alphabetically (ascending)");
    readConfig(config, "removeDuplicates",   removeDuplicates, Config::DEFAULT,  "0", "remove duplicate entries (order is preserved)");
    readConfig(config, "startIndex",         startIndex,       Config::DEFAULT,  "0", "start at element startIndex (counting from 0)");
    readConfig(config, "count",              countElements,    Config::OPTIONAL, "",  "use count elements (default: use all)");
    if(sort)
    {
      logWarningOnce<<"DEPRECATED since 2025-09-27: In Loop->FileStringList->sort . Please use SortAndRemoveDuplicates instead."<<Log::endl;
      std::sort(strings.begin(), strings.end());
    }
    if(removeDuplicates)
    {
      logWarningOnce<<"DEPRECATED since 2025-09-27: In Loop->FileStringList->removeDuplicates. Please use SortAndRemoveDuplicates instead."<<Log::endl;
      std::unordered_set<std::string> set;
      auto end = std::copy_if(strings.begin(), strings.end(), strings.begin(), [&set](auto &s) {return set.insert(s).second;});
      strings.erase(end, strings.end());
    }
    if((startIndex != 0) || (countElements != MAX_UINT))
      logWarningOnce<<"DEPRECATED since 2025-09-27: In Loop->FileStringList->startIndex and ->count. Please use a condition instead."<<Log::endl;
    strings.erase(strings.begin(), strings.begin()+std::min(startIndex, strings.size()));
    strings.erase(strings.begin()+std::min(countElements, strings.size()), strings.end());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Bool LoopFileStringList::iteration(VariableList &varList)
{
  if(index() >= count())
    return FALSE;

  if(!nameString.empty()) varList.setVariable(nameString, strings.at(index()));
  if(!nameIndex.empty())  varList.setVariable(nameIndex,  index());
  if(!nameCount.empty())  varList.setVariable(nameCount,  count());

  return checkCondition(varList);
}

/***********************************************/

#endif
