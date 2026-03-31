/***********************************************/
/**
* @file loopManualTable.h
*
* @brief DEPRECATED since 2025-09-27. Use LoopStringTable instead.
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
DEPRECATED since 2025-09-27. Use LoopStringTable instead.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "inputOutput/logging.h"
#include "classes/loop/loop.h"

/***** CLASS ***********************************/

/** @brief DEPRECATED since 2025-09-27. Use LoopStringTable instead.
* @ingroup loopGroup
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
  std::vector<std::string> nameString;
  std::string              nameIndex, nameCount;

public:
  LoopManualTable(Config &config);

  UInt count() const override {return rows.size();}
  Bool iteration(VariableList &varList) override;
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

inline LoopManualTable::LoopManualTable(Config &config)
{
  try
  {
    std::vector<CellVector> columns;

    std::string choice;
    if(readConfigChoice(config, "table", choice, Config::MUSTSET, "", "define table by rows/columns"))
    {
      if(readConfigChoiceElement(config, "rowWise",    choice, "define table by rows"))
        readConfig(config, "row",    rows,    Config::MUSTSET, "", "define table by rows");
      if(readConfigChoiceElement(config, "columnWise", choice, "define table by columns"))
        readConfig(config, "column", columns, Config::MUSTSET, "", "define table by columns");
      endChoice(config);
    }
    readConfig(config, "variableLoopString", nameString, Config::MUSTSET,  "loopString", "1. variable name for the 1. column, next variable name for the 2. column, ... ");
    readConfig(config, "variableLoopIndex",  nameIndex,  Config::OPTIONAL, "",           "variable with index of current iteration (starts with zero)");
    readConfig(config, "variableLoopCount",  nameCount,  Config::OPTIONAL, "",           "variable with total number of iterations");
    readConfigCondition(config);
    if(isCreateSchema(config)) return;

    logWarningOnce<<"DEPRECATED since 2025-09-27: Loop->ManualTable. Use Loop->StringTable instead."<<Log::endl;

    // transpose columns
    for(UInt i=0; i<columns.size(); i++)
      for(UInt k=0; k<columns.at(i).cells.size(); k++)
      {
        rows.resize(k+1);
        rows.at(k).cells.push_back(columns.at(i).cells.at(k));
      }

    if(rows.size() > 1)
    {
      if(nameString.size() > rows.at(0).cells.size())
        throw(Exception("More variables than columns in the table"));
      for(UInt i=1; i<rows.size(); i++)
        if(rows.at(i).cells.size() != rows.at(0).cells.size())
          throw(Exception("Varying number of columns in the table"));
     }
   }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Bool LoopManualTable::iteration(VariableList &varList)
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
