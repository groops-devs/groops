/***********************************************/
/**
* @file loopFileLines.h
*
* @brief Loop over lines of a text file.
*
* @author Torsten Mayer-Guerr
* @date 2023-05-31
*
*/
/***********************************************/

#ifndef __GROOPS_LOOPFILELINES__
#define __GROOPS_LOOPFILELINES__

// Latex documentation
#ifdef DOCSTRING_Loop
static const char *docstringLoopFileLines = R"(
\subsection{FileLines}
Loop over lines of a text file.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "inputOutput/file.h"
#include "classes/loop/loop.h"
#include <unordered_set>

/***** CLASS ***********************************/

/** @brief Loop over lines of a text file.
* @ingroup loopGroup
* @see Loop */
class LoopFileLines : public Loop
{
  std::string nameLine, nameIndex, nameCount;
  std::vector<std::string> lines;

public:
  LoopFileLines(Config &config);

  UInt count() const override {return lines.size();}
  Bool iteration(VariableList &varList) override;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline LoopFileLines::LoopFileLines(Config &config)
{
  try
  {
    std::vector<FileName> fileNames;
    UInt startIndex;
    UInt countElements = MAX_UINT;
    Bool sort;
    Bool removeDuplicates;

    readConfig(config, "inputfile",         fileNames,        Config::MUSTSET,  "",         "simple ASCII file with lines");
    readConfig(config, "sort",              sort,             Config::DEFAULT,  "0",        "sort lines alphabetically (ascending)");
    readConfig(config, "removeDuplicates",  removeDuplicates, Config::DEFAULT,  "0",        "remove duplicate lines (order is preserved)");
    readConfig(config, "startIndex",        startIndex,       Config::DEFAULT,  "0",        "start at element startIndex (counting from 0)");
    readConfig(config, "count",             countElements,    Config::OPTIONAL, "",         "use number of lines only (default: use all)");
    readConfig(config, "variableLoopLine",  nameLine,         Config::OPTIONAL, "loopLine", "name of the variable to be replaced");
    readConfig(config, "variableLoopIndex", nameIndex,        Config::OPTIONAL, "",         "variable with index of current iteration (starts with zero)");
    readConfig(config, "variableLoopCount", nameCount,        Config::OPTIONAL, "",         "variable with total number of iterations");
    readConfigCondition(config);
    if(isCreateSchema(config)) return;

    for(auto &fileName : fileNames)
    {
      InFile file(fileName);
      std::string line;
      while(std::getline(file, line))
        lines.push_back(line);
    }

    if(sort)
      std::sort(lines.begin(), lines.end());

    if(removeDuplicates)
    {
      std::unordered_set<std::string> set;
      lines.erase(std::copy_if(lines.begin(), lines.end(), lines.begin(), [&set](auto &s) {return set.insert(s).second;}), lines.end());
    }

    lines.erase(lines.begin(), lines.begin()+std::min(startIndex, lines.size()));
    lines.erase(lines.begin()+std::min(countElements, lines.size()), lines.end());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Bool LoopFileLines::iteration(VariableList &varList)
{
  if(index() >= count())
    return FALSE;

  if(!nameLine.empty())  varList.setVariable(nameLine,  lines.at(index()));
  if(!nameIndex.empty()) varList.setVariable(nameIndex, index());
  if(!nameCount.empty()) varList.setVariable(nameCount, count());

  return checkCondition(varList);
}

/***********************************************/

#endif
