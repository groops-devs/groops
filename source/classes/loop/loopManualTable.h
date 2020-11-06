/***********************************************/
/**
* @file loopManualTable.h
*
* @brief Loop over rows of a table containing strings.
*
* @author Sebastian Strasser
* @date 2017-02-20
*
*/
/***********************************************/

#ifndef __GROOPS_LOOPMANUALTABLE__
#define __GROOPS_LOOPMANUALTABLE__

// Latex documentation
#ifdef DOCSTRING_Loop
static const char *docstringLoopManualTable = R"(
\subsection{ManualTable}
Loop over rows of a table containing strings.
The table can be defined \config{rowWise} or \config{columnWise}.
Each row/column must have the same number of cells.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "inputOutput/logging.h"
#include "classes/loop/loop.h"

/***** CLASS ***********************************/

/** @brief Loop over rows of a table containing strings.
* @ingroup LoopGroup
* @see Loop */
class LoopManualTable : public Loop
{
public:
  class CellVector
  {
  public:
    std::vector<std::string> cells;
  };

private:
  std::vector<CellVector>  rows;
  std::vector<CellVector>  columns;
  std::vector<std::string> nameString;
  std::string              nameIndex, nameCount;

public:
  LoopManualTable(Config &config);

  UInt count() const override;
  void init(VariableList &varList) override;
  void setValues(VariableList &varList) override;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

template<> Bool readConfig(Config &config, const std::string &name, LoopManualTable::CellVector &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  if(!readConfigSequence(config, name, mustSet, defaultValue, annotation))
    return FALSE;
  readConfig(config, "cell", var.cells, Config::MUSTSET, "", "explicit list of cells in row/column");
  endSequence(config);
  return TRUE;
}

/***********************************************/

inline LoopManualTable::LoopManualTable(Config &config) : Loop()
{
  try
  {
    std::string choice;
    readConfigChoice(config, "table",  choice, Config::MUSTSET, "", "define table by rows/columns");
    if(readConfigChoiceElement(config, "rowWise",    choice, "define table by rows"))
      readConfig(config, "row",              rows,       Config::MUSTSET,   "",           "define table by rows");
    if(readConfigChoiceElement(config, "columnWise", choice, "define table by columns"))
      readConfig(config, "column",           columns,    Config::MUSTSET,   "",           "define table by columns");
    endChoice(config);
    readConfig(config, "variableLoopString", nameString, Config::MUSTSET,   "loopString", "1. variable name for the 1. column, next variable name for the 2. column, ... ");
    readConfig(config, "variableLoopIndex",  nameIndex,  Config::OPTIONAL,  "",           "variable with index of current iteration (starts with zero)");
    readConfig(config, "variableLoopCount",  nameCount,  Config::OPTIONAL,  "",           "variable with total number of iterations");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline UInt LoopManualTable::count() const
{
  return (columns.size() ? columns.at(0).cells.size() : rows.size());
}

/***********************************************/

inline void LoopManualTable::init(VariableList &varList)
{
  try
  {
    if(rows.size() > 1 || columns.size() > 1)
    {
      const UInt rowCount    = (columns.size() ? columns.at(0).cells.size() : rows.size());
      const UInt columnCount = (rows.size()    ? rows.at(0).cells.size()    : columns.size());
      if(nameString.size()>columnCount)
        throw(Exception("More variables than columns in the table"));
      for(UInt i = 1; i < rows.size(); i++)
        if(rows.at(i).cells.size() != columnCount)
          throw(Exception("Varying number of columns in the table"));
      for(UInt i = 1; i < columns.size(); i++)
        if(columns.at(i).cells.size() != rowCount)
          throw(Exception("Varying number of rows in the table"));
    }

    if(index == NULLINDEX)
    {
      for(UInt i=0; i<nameString.size(); i++)
        if(!nameString.at(i).empty())
          addVariable(nameString.at(i), varList);
      if(!nameIndex.empty()) addVariable(nameIndex, varList);
      if(!nameCount.empty()) addVariable(nameCount, varList);
    }
    index = 0;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void LoopManualTable::setValues(VariableList &varList)
{
  for(UInt i=0; i<nameString.size(); i++)
    if(!nameString.at(i).empty())
    {
      if(rows.size())
        varList[nameString.at(i)]->setValue(rows.at(index).cells.at(i));
      else if(columns.size())
        varList[nameString.at(i)]->setValue(columns.at(i).cells.at(index));
    }
  if(!nameIndex.empty())  varList[nameIndex]->setValue(index);
  if(!nameCount.empty())  varList[nameCount]->setValue(count());
}


/***********************************************/

#endif
