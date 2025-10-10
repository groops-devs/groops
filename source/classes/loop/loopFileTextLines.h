/***********************************************/
/**
* @file loopFileTextLines.h
*
* @brief Loop over lines of a text file.
*
* @author Torsten Mayer-Guerr
* @date 2023-05-31
*
*/
/***********************************************/

#ifndef __GROOPS_LOOPFILETEXTLINES__
#define __GROOPS_LOOPFILETEXTLINES__

// Latex documentation
#ifdef DOCSTRING_Loop
static const char *docstringLoopFileTextLines = R"(
\subsection{FileTextLines}\label{loopType:fileTextLines}
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
class LoopFileTextLines : public Loop
{
  std::vector<std::string> lines;
  std::string nameLine, nameIndex, nameCount;

public:
  LoopFileTextLines(Config &config);

  UInt count() const override {return lines.size();}
  Bool iteration(VariableList &varList) override;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline LoopFileTextLines::LoopFileTextLines(Config &config)
{
  try
  {
    std::vector<FileName> fileNames;
    UInt startIndex, countElements = MAX_UINT;

    readConfig(config, "inputfile",         fileNames,     Config::MUSTSET,  "",         "simple text file with lines");
    readConfig(config, "startIndex",        startIndex,    Config::DEFAULT,  "0",        "start at element startIndex (counting from 0)");
    readConfig(config, "count",             countElements, Config::OPTIONAL, "",         "use number of loops only (default: use all)");
    readConfig(config, "variableLoopLine",  nameLine,      Config::OPTIONAL, "loopLine", "name of the variable to be replaced");
    readConfig(config, "variableLoopIndex", nameIndex,     Config::OPTIONAL, "",         "variable with index of current iteration (starts with zero)");
    readConfig(config, "variableLoopCount", nameCount,     Config::OPTIONAL, "",         "variable with total number of iterations");
    readConfigCondition(config);
    if(isCreateSchema(config)) return;

    for(auto &fileName : fileNames)
    {
      InFile file(fileName);
      std::string line;
      while(std::getline(file, line))
        lines.push_back(line);
    }

    lines.erase(lines.begin(), lines.begin()+std::min(startIndex, lines.size()));
    lines.erase(lines.begin()+std::min(countElements, lines.size()), lines.end());

    // DEPRECTATED
    // -----------
    // read deprecated config elements (as after if(isCreateSchema(config)) not shown in the GUI and docu)
    Bool sort;
    Bool removeDuplicates;
    readConfig(config, "sort",             sort,             Config::DEFAULT,  "0", "sort lines alphabetically (ascending)");
    readConfig(config, "removeDuplicates", removeDuplicates, Config::DEFAULT,  "0", "remove duplicate lines (order is preserved)");
    if(sort)
    {
      logWarningOnce<<"DEPRECATED since 2025-09-27: Loop->FileTextLines->sort. Please use Loop->SortAndRemoveDuplicates instead."<<Log::endl;
      std::sort(lines.begin(), lines.end());
    }
    if(removeDuplicates)
    {
      logWarningOnce<<"DEPRECATED since 2025-09-27: Loop->FileTextLines->removeDuplicates. Please use Loop->SortAndRemoveDuplicates instead."<<Log::endl;
      std::unordered_set<std::string> set;
      lines.erase(std::copy_if(lines.begin(), lines.end(), lines.begin(), [&set](auto &s) {return set.insert(s).second;}), lines.end());
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Bool LoopFileTextLines::iteration(VariableList &varList)
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
