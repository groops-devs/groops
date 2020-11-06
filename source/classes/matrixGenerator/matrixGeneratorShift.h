/***********************************************/
/**
* @file matrixGeneratorShift.h
*
* @brief Shift start row and start column of a matrix.
*
* @author Torsten Mayer-Guerr
* @date 2017-09-01
*
*/
/***********************************************/

#ifndef __GROOPS_MATRIXCGENERATORSHIFT__
#define __GROOPS_MATRIXCGENERATORSHIFT__

// Latex documentation
#ifdef DOCSTRING_MatrixGenerator
static const char *docstringMatrixGeneratorShift = R"(
\subsection{Shift}
Shift start row and start column of a matrix.
In other words: zero lines and columns are inserted at the beginning of the matrix.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "matrixGenerator.h"

/***** CLASS ***********************************/

/** @brief Shift of a matrix.
* @ingroup matrixGeneratorGroup
* @see MatrixGenerator */
class MatrixGeneratorShift : public MatrixGeneratorBase
{
  MatrixGeneratorPtr    matrix;
  ExpressionVariablePtr exprStartRow, exprStartCol;

public:
  MatrixGeneratorShift(Config &config);
  void compute(Matrix &A, UInt &startRow, UInt &startCol);
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline MatrixGeneratorShift::MatrixGeneratorShift(Config &config) : MatrixGeneratorBase(config)
{
  try
  {
    readConfig(config, "matrix",      matrix,       Config::MUSTSET,  "",  "");
    readConfig(config, "startRow",    exprStartRow, Config::DEFAULT,  "0", "start row (variables: rowsBefore, columnsBefore, rows, columns)");
    readConfig(config, "startColumn", exprStartCol, Config::DEFAULT,  "0", "start column (variables: rowsBefore, columnsBefore, rows, columns)");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void MatrixGeneratorShift::compute(Matrix &A, UInt &startRow, UInt &startCol)
{
  try
  {
    A = matrix->compute();

    addVariable("rows",    static_cast<Double>(A.rows()),    varList);
    addVariable("columns", static_cast<Double>(A.columns()), varList);
    startRow = static_cast<UInt>(exprStartRow->evaluate(varList));
    startCol = static_cast<UInt>(exprStartCol->evaluate(varList));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
