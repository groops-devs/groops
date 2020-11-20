/***********************************************/
/**
* @file loopFileAscii.h
*
* @brief Loop over list of strings from files.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2017-01-27
*
*/
/***********************************************/

#ifndef __GROOPS_LOOPFILEASCII__
#define __GROOPS_LOOPFILEASCII__

// Latex documentation
#ifdef DOCSTRING_Loop
static const char *docstringLoopFileAscii = R"(
\subsection{FileAscii}
Loop over list of strings from files.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "files/fileStringTable.h"
#include "classes/loop/loop.h"
#include <unordered_set>

/***** CLASS ***********************************/

/** @brief Loop over list of strings from files.
* @ingroup LoopGroup
* @see Loop */
class LoopFileAscii : public Loop
{
  std::string nameString, nameIndex, nameCount;
  std::vector<std::string> strings;
  UInt index;

public:
  LoopFileAscii(Config &config);

  UInt count() const override {return strings.size();}
  Bool iteration(VariableList &varList) override;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline LoopFileAscii::LoopFileAscii(Config &config)
{
  try
  {
    std::vector<FileName> fileName;
    UInt startIndex;
    UInt countElements;
    Bool sort;
    Bool removeDuplicates;

    readConfig(config, "inputfile",          fileName,         Config::MUSTSET,   "",           "simple ASCII file with strings (separated by whitespace)");
    readConfig(config, "sort",               sort,             Config::DEFAULT,   "0",          "sort entries alphabetically (ascending)");
    readConfig(config, "removeDuplicates",   removeDuplicates, Config::DEFAULT,   "0",          "remove duplicate entries (order is preserved)");
    readConfig(config, "startIndex",         startIndex,       Config::DEFAULT,   "0",          "start at element startIndex (counting from 0)");
    readConfig(config, "count",              countElements,    Config::OPTIONAL, "",           "use count elements (default: use all)");
    readConfig(config, "variableLoopString", nameString,       Config::OPTIONAL, "loopString", "name of the variable to be replaced");
    readConfig(config, "variableLoopIndex",  nameIndex,        Config::OPTIONAL, "",           "variable with index of current iteration (starts with zero)");
    readConfig(config, "variableLoopCount",  nameCount,        Config::OPTIONAL, "",           "variable with total number of iterations");
    if(isCreateSchema(config)) return;

    for(UInt i=0; i<fileName.size(); i++)
    {
      std::vector<std::string> str;
      readFileStringList(fileName.at(i), str);
      strings.insert(strings.end(), str.begin(), str.end());
    }

    if(sort)
      std::sort(strings.begin(), strings.end());

    if(removeDuplicates)
    {
      std::unordered_set<std::string> set;
      auto end = std::copy_if(strings.begin(), strings.end(), strings.begin(), [&set](auto &s) {return set.insert(s).second;});
      strings.erase(end, strings.end());
    }

    strings.erase(strings.begin(), strings.begin()+std::min(startIndex, strings.size()));
    strings.erase(strings.begin()+std::min(countElements, strings.size()), strings.end());

    index = 0;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Bool LoopFileAscii::iteration(VariableList &varList)
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
