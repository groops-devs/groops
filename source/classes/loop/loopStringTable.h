/***********************************************/
/**
* @file loopStringTable.h
*
* @brief Loop over rows of a table containing strings.
*
* @author Sebastian Strasser
* @date 2017-02-20
*
*/
/***********************************************/

#ifndef __GROOPS_LOOPSTRINGTABLE__
#define __GROOPS_LOOPSTRINGTABLE__

// Latex documentation
#ifdef DOCSTRING_Loop
static const char *docstringLoopStringTable = R"(
\subsection{StringTable}\label{loopType:stringTable}
Loop over \config{row}s of a table containing strings.
Each row must have the same number of cells.
For each column an extra \config{variableLoopString} can be defined.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "inputOutput/logging.h"
#include "classes/loop/loop.h"

/***** CLASS ***********************************/

/** @brief Loop over rows of a table containing strings.
* @ingroup loopGroup
* @see Loop */
class LoopStringTable : public Loop
{
public:
  class CellVector
  {
  public:
    std::vector<std::string> cells;
  };

private:
  std::vector<CellVector>  rows;
  std::vector<std::string> nameString;
  std::string              nameIndex, nameCount;

public:
  LoopStringTable(Config &config);

  UInt count() const override {return rows.size();}
  Bool iteration(VariableList &varList) override;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

template<> Bool readConfig(Config &config, const std::string &name, LoopStringTable::CellVector &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  if(!readConfigSequence(config, name, mustSet, defaultValue, annotation))
    return FALSE;
  readConfig(config, "cell", var.cells, Config::MUSTSET, "", "list of columns in a row");
  endSequence(config);
  return TRUE;
}

/***********************************************/

inline LoopStringTable::LoopStringTable(Config &config)
{
  try
  {
    Bool transpose;

    readConfig(config, "row",                rows,       Config::OPTIONAL, "",           "rows of a table");
    readConfig(config, "transpose",          transpose,  Config::DEFAULT,  "0",          "loop over columns instead of rows");
    readConfig(config, "variableLoopString", nameString, Config::MUSTSET,  "loopString", "1. variable name for the 1. column, next variable name for the 2. column, ... ");
    readConfig(config, "variableLoopIndex",  nameIndex,  Config::OPTIONAL, "",           "variable with index of current iteration (starts with zero)");
    readConfig(config, "variableLoopCount",  nameCount,  Config::OPTIONAL, "",           "variable with total number of iterations");
    readConfigCondition(config);
    if(isCreateSchema(config)) return;

    // transpose table
    if(rows.size() && transpose)
    {
      const UInt cols = std::max_element(rows.begin(), rows.end(), [](auto &s1, auto &s2){return s1.cells.size() < s2.cells.size();})->cells.size();
      std::vector<CellVector> tmp(cols);
      for(UInt i=0; i<rows.size(); i++)
        for(UInt k=0; k<rows.at(i).cells.size(); k++)
          tmp.at(k).cells.push_back(rows.at(i).cells.at(k));
      rows = std::move(tmp);
    }

    if(rows.size() && (nameString.size() > std::min_element(rows.begin(), rows.end(), [](auto &s1, auto &s2){return s1.cells.size() < s2.cells.size();})->cells.size()))
      throw(Exception("More variables than columns in the table"));
   }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Bool LoopStringTable::iteration(VariableList &varList)
{
  if(index() >= count())
    return FALSE;

  for(UInt i=0; i<nameString.size(); i++)
    if(!nameString.at(i).empty())
      varList.setVariable(nameString.at(i), rows.at(index()).cells.at(i));
  if(!nameIndex.empty()) varList.setVariable(nameIndex, index());
  if(!nameCount.empty()) varList.setVariable(nameCount, count());

  return checkCondition(varList);
}

/***********************************************/

#endif
