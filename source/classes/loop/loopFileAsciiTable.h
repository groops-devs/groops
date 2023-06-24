/***********************************************/
/**
* @file loopFileAsciiTable.h
*
* @brief Loop over rows of a table containing strings.
*
* @author Norbert Zehentner
* @author Sebastian Strasser
* @date 2017-02-17
*
*/
/***********************************************/

#ifndef __GROOPS_LOOPFILEASCIITABLE__
#define __GROOPS_LOOPFILEASCIITABLE__

// Latex documentation
#ifdef DOCSTRING_Loop
static const char *docstringLoopFileAsciiTable = R"(
\subsection{FileAsciiTable}
Loop over rows of a table containing strings.
Each row must have the same number of columns.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "files/fileStringTable.h"
#include "classes/loop/loop.h"

/***** CLASS ***********************************/

/** @brief Loop over rows of a table containing strings.
* @ingroup LoopGroup
* @see Loop */
class LoopFileAsciiTable : public Loop
{
  std::vector<std::vector<std::string>>  stringTable;
  std::vector<std::string>               nameString;
  std::string                            nameIndex, nameCount;
  UInt                                   index;

public:
  LoopFileAsciiTable(Config &config);

  UInt count() const override {return stringTable.size();}
  Bool iteration(VariableList &varList) override;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline LoopFileAsciiTable::LoopFileAsciiTable(Config &config)
{
  try
  {
    std::vector<FileName> fileNames;
    UInt     startRow, countRows = MAX_UINT;

    readConfig(config, "inputfile",          fileNames,  Config::MUSTSET,  "",           "simple ASCII file with multiple columns (separated by whitespace)");
    readConfig(config, "startLine",          startRow,   Config::DEFAULT,  "0",          "start at line startLine (counting from 0)");
    readConfig(config, "countLines",         countRows,  Config::OPTIONAL, "",           "read count lines (default: all)");
    readConfig(config, "variableLoopString", nameString, Config::MUSTSET,  "loopString", "1. variable name for the 1. column, next variable name for the 2. column, ... ");
    readConfig(config, "variableLoopIndex",  nameIndex,  Config::OPTIONAL, "",           "variable with index of current iteration (starts with zero)");
    readConfig(config, "variableLoopCount",  nameCount,  Config::OPTIONAL, "",           "variable with total number of iterations");
    if(isCreateSchema(config)) return;

    for(auto &fileName : fileNames)
      readFileStringTable(fileName, stringTable);

    stringTable.erase(stringTable.begin(), stringTable.begin()+std::min(startRow, stringTable.size()));
    stringTable.erase(stringTable.begin()+std::min(countRows, stringTable.size()), stringTable.end());

    if(stringTable.size()>1)
    {
      const UInt columns = stringTable.at(0).size();
      if(nameString.size() > columns)
        throw(Exception("More variables than columns in the table"));
      for(UInt i=1; i<stringTable.size(); i++)
        if(stringTable.at(i).size() != columns)
          throw(Exception("Varying number of columns in input file(s)"));
    }

    index = 0;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Bool LoopFileAsciiTable::iteration(VariableList &varList)
{
  if(index >= count())
    return FALSE;

  for(UInt i=0; i<nameString.size(); i++)
    if(!nameString.at(i).empty())
      varList.setVariable(nameString.at(i), stringTable.at(index).at(i));
  if(!nameIndex.empty()) varList.setVariable(nameIndex, index);
  if(!nameCount.empty()) varList.setVariable(nameCount, count());

  index++;
  return TRUE;
}

/***********************************************/

#endif
