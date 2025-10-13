/***********************************************/
/**
* @file loopFileStringTable.h
*
* @brief Loop over rows of a table containing strings.
*
* @author Norbert Zehentner
* @author Sebastian Strasser
* @date 2017-02-17
*
*/
/***********************************************/

#ifndef __GROOPS_LOOPFILESTRINGTABLE__
#define __GROOPS_LOOPFILESTRINGTABLE__

// Latex documentation
#ifdef DOCSTRING_Loop
static const char *docstringLoopFileStringTable = R"(
\subsection{FileStringTable}\label{loopType:fileStringTable}
Loop over rows of a \file{table}{stringTable} containing strings.
Each row must have the same number of columns.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "files/fileStringTable.h"
#include "classes/loop/loop.h"

/***** CLASS ***********************************/

/** @brief Loop over rows of a table containing strings.
* @ingroup loopGroup
* @see Loop */
class LoopFileStringTable : public Loop
{
  std::vector<std::vector<std::string>>  stringTable;
  std::vector<std::string>               nameString;
  std::string                            nameIndex, nameCount;

public:
  LoopFileStringTable(Config &config);

  UInt count() const override {return stringTable.size();}
  Bool iteration(VariableList &varList) override;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline LoopFileStringTable::LoopFileStringTable(Config &config)
{
  try
  {
    std::vector<FileName> fileNames;
    Bool transpose;

    readConfig(config, "inputfile",          fileNames,  Config::MUSTSET,  "",           "string table file with multiple columns");
    readConfig(config, "transpose",          transpose,  Config::DEFAULT,  "0",          "loop over columns instead of rows");
    readConfig(config, "variableLoopString", nameString, Config::MUSTSET,  "loopString", "1. variable name for the 1. column, next variable name for the 2. column, ... ");
    readConfig(config, "variableLoopIndex",  nameIndex,  Config::OPTIONAL, "",           "variable with index of current iteration (starts with zero)");
    readConfig(config, "variableLoopCount",  nameCount,  Config::OPTIONAL, "",           "variable with total number of iterations");
    readConfigCondition(config);
    if(isCreateSchema(config)) return;

    for(auto &fileName : fileNames)
      readFileStringTable(fileName, stringTable);

    // transpose table
    if(stringTable.size() && transpose)
    {
      const UInt cols = std::max_element(stringTable.begin(), stringTable.end(), [](auto &s1, auto &s2){return s1.size() < s2.size();})->size();
      std::vector<std::vector<std::string>> tmp(cols);
      for(UInt i=0; i<stringTable.size(); i++)
        for(UInt k=0; k<stringTable.at(i).size(); k++)
          tmp.at(k).push_back(stringTable.at(i).at(k));
      stringTable = std::move(tmp);
    }

    if(stringTable.size() && (nameString.size() > std::min_element(stringTable.begin(), stringTable.end(), [](auto &s1, auto &s2){return s1.size() < s2.size();})->size()))
      throw(Exception("More variables than minimum columns in the table"));

    // DEPRECTATED
    // -----------
    // read deprecated config elements (as after if(isCreateSchema(config)) not shown in the GUI and docu)
    UInt startRow, countRows = MAX_UINT;
    readConfig(config, "startLine",  startRow,   Config::DEFAULT,  "0", "start at line startLine (counting from 0)");
    readConfig(config, "countLines", countRows,  Config::OPTIONAL, "",  "read count lines (default: all)");
    if((startRow != 0) || (countRows != MAX_UINT))
      logWarningOnce<<"DEPRECATED since 2025-09-27: In Loop->FileStringTable->startIndex and ->count. Please use a condition instead."<<Log::endl;
    stringTable.erase(stringTable.begin(), stringTable.begin()+std::min(startRow, stringTable.size()));
    stringTable.erase(stringTable.begin()+std::min(countRows, stringTable.size()), stringTable.end());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Bool LoopFileStringTable::iteration(VariableList &varList)
{
  if(index() >= count())
    return FALSE;

  for(UInt i=0; i<nameString.size(); i++)
    if(!nameString.at(i).empty())
      varList.setVariable(nameString.at(i), stringTable.at(index()).at(i));
  if(!nameIndex.empty()) varList.setVariable(nameIndex, index());
  if(!nameCount.empty()) varList.setVariable(nameCount, count());

  return checkCondition(varList);
}

/***********************************************/

#endif
