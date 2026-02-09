/***********************************************/
/**
* @file loopSortAndRemoveDuplicates.h
*
* @brief Sort and removes duplicates from a loop.
*
* @author Torsten Mayer-Guerr
* @date 2025-09-25
*
*/
/***********************************************/

#ifndef __GROOPS_LOOPSORTANDREMOVEDUPLICATES__
#define __GROOPS_LOOPSORTANDREMOVEDUPLICATES__

// Latex documentation
#ifdef DOCSTRING_Loop
static const char *docstringLoopSortAndRemoveDuplicates = R"(
\subsection{SortAndRemoveDuplicates}\label{loopType:sortAndRemoveDuplicates}
Perform the \configClass{loop}{loopType} in the alphabetically
order defined by the evaluated \config{sortString} for each loop step.
So the string must contain loop variables. If \config{sortString}
is empty, no sorting will take place.

Example: The \config{sortString}=\verb|{loopTime:%m}| of a time series
sorts the times in ascending order by month.

The same principle is used to remove duplicates. If different loop steps
evaluates \config{removeDuplicatesString} to the same string,
only the first loop step is executed.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "parser/stringParser.h"
#include "classes/loop/loop.h"
#include <unordered_set>

/***** CLASS ***********************************/

/** @brief Loop over lines of a text file.
* @ingroup loopGroup
* @see Loop */
class LoopSortAndRemoveDuplicates : public Loop
{
  class Loop
  {
  public:
    VariableList varList;
    std::string sortString, duplicatesString;
  };

  std::vector<Loop> loops;
  std::string nameIndex, nameCount;

public:
  LoopSortAndRemoveDuplicates(Config &config);

  UInt count() const override {return loops.size();}
  Bool iteration(VariableList &varList) override;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline LoopSortAndRemoveDuplicates::LoopSortAndRemoveDuplicates(Config &config)
{
  try
  {
    LoopPtr loopPtr;
    std::pair<std::string, VariableList> sortString, duplicatesString;
    Bool    descending;

    readConfig(config, "loop",                   loopPtr,          Config::MUSTSET,  "",  "");
    readConfig(config, "sortString",             sortString,       Config::OPTIONAL, "",  "use {loopVariables}, sort alphabetically");
    readConfig(config, "descending",             descending,       Config::DEFAULT,  "0", "sorting descending instead of ascending");
    readConfig(config, "removeDuplicatesString", duplicatesString, Config::OPTIONAL, "",  "use {loopVariables}, remove duplicates (order is preserved)");
    readConfig(config, "variableLoopIndex",      nameIndex,        Config::OPTIONAL, "",  "variable with index of current iteration (starts with zero)");
    readConfig(config, "variableLoopCount",      nameCount,        Config::OPTIONAL, "",  "variable with total number of iterations");
    readConfigCondition(config);
    if(isCreateSchema(config)) return;

    // perform the loops
    // -----------------
    VariableList varList;
    for(;;)
    {
      loops.resize(loops.size()+1);
      if(!loopPtr->iteration(varList))
        break;
      loops.back().varList = varList;
      // evaluate sort string & duplicate string
      loops.back().sortString       = StringParser::parse(sortString.first,       VariableList(sortString.second)       += loops.back().varList);
      loops.back().duplicatesString = StringParser::parse(duplicatesString.first, VariableList(duplicatesString.second) += loops.back().varList);
    }
    loops.pop_back(); // remove last unused

    // sort
    // ----
    if(!sortString.first.empty())
    {
      if(descending)
        std::stable_sort(loops.begin(), loops.end(), [](auto &l1, auto &l2){return l1.sortString > l2.sortString;});
      else
        std::stable_sort(loops.begin(), loops.end(), [](auto &l1, auto &l2){return l1.sortString < l2.sortString;});
    }

    // remove duplicates
    // -----------------
    if(!duplicatesString.first.empty())
    {
      std::unordered_set<std::string> set;
      loops.erase(std::copy_if(loops.begin(), loops.end(), loops.begin(), [&set](auto &loop) {return set.insert(loop.duplicatesString).second;}), loops.end());
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Bool LoopSortAndRemoveDuplicates::iteration(VariableList &varList)
{
  if(index() >= count())
    return FALSE;

  varList += loops.at(index()).varList;
  if(!nameIndex.empty()) varList.setVariable(nameIndex, index());
  if(!nameCount.empty()) varList.setVariable(nameCount, count());

  return checkCondition(varList);
}

/***********************************************/

#endif
