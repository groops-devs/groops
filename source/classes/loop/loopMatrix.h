/***********************************************/
/**
* @file loopMatrix.h
*
* @brief Loop over rows of a matrix.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2017-01-27
*
*/
/***********************************************/

#ifndef __GROOPS_LOOPMATRIX__
#define __GROOPS_LOOPMATRIX__

// Latex documentation
#ifdef DOCSTRING_Loop
static const char *docstringLoopMatrix = R"(
\subsection{Matrix}
Loop over rows of a matrix. To define the loop variables the standard
data variables of the matrix are available, see~\reference{dataVariables}{general.parser:dataVariables}.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "parser/dataVariables.h"
#include "files/fileMatrix.h"
#include "classes/loop/loop.h"

/***** CLASS ***********************************/

/** @brief Loop over rows of a matrix.
* @ingroup LoopGroup
* @see Loop */
class LoopMatrix : public Loop
{
  Matrix       A;
  mutable VariableList varListMatrix;
  std::vector<ExpressionVariablePtr> variableExpr;
  std::string nameIndex, nameCount;

  FileName fileName;
  ExpressionVariablePtr variableStartRow;
  ExpressionVariablePtr variableCountRows;
  Bool transpose;

public:
  LoopMatrix(Config &config);

  UInt count() const override;
  void init(VariableList &varList) override;
  void setValues(VariableList &varList) override;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline LoopMatrix::LoopMatrix(Config &config) : Loop()
{
  try
  {
    readConfig(config, "inputfile",          fileName,          Config::MUSTSET,   "",                 "");
    readConfig(config, "transpose",          transpose,         Config::DEFAULT,   "0",                "effectively loop over columns");
    readConfig(config, "startRow",           variableStartRow,  Config::DEFAULT,   "0",                "start at this row (variable: rows)");
    readConfig(config, "countRows",          variableCountRows, Config::DEFAULT,   "rows",             "use this many rows (variable: rows)");
    readConfig(config, "variableLoop",       variableExpr,      Config::MUSTSET,   "loopNumber=data0", "define a variable by name = expression (input columns are named data0, data1, ...)");
    readConfig(config, "variableLoopIndex",  nameIndex,         Config::OPTIONAL, "",                 "variable with index of current iteration (starts with zero)");
    readConfig(config, "variableLoopCount",  nameCount,         Config::OPTIONAL, "",                 "variable with total number of iterations");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline UInt LoopMatrix::count() const
{
  return A.rows();
}

/***********************************************/

inline void LoopMatrix::init(VariableList &varList)
{
  try
  {
    varListMatrix = varList;

    readFileMatrix(fileName(varList), A);
    if(transpose)
      A = A.trans();

    addVariable("rows", A.rows(), varListMatrix);
    A = A.row(static_cast<UInt>(variableStartRow->evaluate(varListMatrix)), static_cast<UInt>(variableCountRows->evaluate(varListMatrix)));

    if(index == NULLINDEX)
    {
      for(UInt i=0; i<variableExpr.size(); i++)
      {
        variableExpr.at(i)->parseVariableName(varListMatrix); // get real variable names, otherwise all named 'variableLoop'
        addVariable(variableExpr.at(i)->name(), varList);
      }
      std::set<std::string> usedVariables;
      for(UInt i=0; i<variableExpr.size(); i++)
        variableExpr.at(i)->usedVariables(varListMatrix, usedVariables);
      addDataVariables(A, varListMatrix, usedVariables);
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

inline void LoopMatrix::setValues(VariableList &varList)
{
  evaluateDataVariables(A, index, varListMatrix);
  for(UInt i=0; i<variableExpr.size(); i++)
    varList[variableExpr.at(i)->name()]->setValue(variableExpr.at(i)->evaluate(varListMatrix));
  if(!nameIndex.empty())  varList[nameIndex]->setValue(index);
  if(!nameCount.empty())  varList[nameCount]->setValue(count());
}

/***********************************************/

#endif
