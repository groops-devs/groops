/***********************************************/
/**
* @file matrixGeneratorSlice.h
*
* @brief Slice of a matrix.
*
* @author Torsten Mayer-Guerr
* @date 2017-09-01
*
*/
/***********************************************/

#ifndef __GROOPS_MATRIXGENERATORSLICE__
#define __GROOPS_MATRIXGENERATORSLICE__

// Latex documentation
#ifdef DOCSTRING_MatrixGenerator
static const char *docstringMatrixGeneratorSlice = R"(
\subsection{Slice}
Slice of a matrix.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "matrixGenerator.h"

/***** CLASS ***********************************/

/** @brief Slice of a matrix.
* @ingroup matrixGeneratorGroup
* @see MatrixGenerator */
class MatrixGeneratorSlice : public MatrixGeneratorBase
{
  MatrixGeneratorPtr    matrix;
  ExpressionVariablePtr exprRow, exprCol;
  ExpressionVariablePtr exprRows, exprCols;

public:
  MatrixGeneratorSlice(Config &config);
  void compute(Matrix &A, UInt &startRow, UInt &startCol);
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline MatrixGeneratorSlice::MatrixGeneratorSlice(Config &config) : MatrixGeneratorBase(config)
{
  try
  {
    readConfig(config, "matrix",      matrix,   Config::MUSTSET,  "",  "");
    readConfig(config, "startRow",    exprRow,  Config::DEFAULT,  "0", "start row of matrix (variables: rowsBefore, columnsBefore, rows, columns)");
    readConfig(config, "startColumn", exprCol,  Config::DEFAULT,  "0", "start column of matrix (variables: rowsBefore, columnsBefore, rows, columns)");
    readConfig(config, "rows",        exprRows, Config::DEFAULT,  "0", "0: until end (variables: rowsBefore, columnsBefore, rows, columns)");
    readConfig(config, "columns",     exprCols, Config::DEFAULT,  "0", "0: until end (variables: rowsBefore, columnsBefore, rows, columns)");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void MatrixGeneratorSlice::compute(Matrix &A, UInt &/*startRow*/, UInt &/*startCol*/)
{
  try
  {
    A = matrix->compute();

    addVariable("rows",    static_cast<Double>(A.rows()),    varList);
    addVariable("columns", static_cast<Double>(A.columns()), varList);
    const UInt row  = static_cast<UInt>(exprRow->evaluate(varList));
    const UInt col  = static_cast<UInt>(exprCol->evaluate(varList));
    const UInt rows = static_cast<UInt>(exprRows->evaluate(varList));
    const UInt cols = static_cast<UInt>(exprCols->evaluate(varList));

    A = A.slice(row, col, rows > 0 ? rows : A.rows()-row, cols > 0 ? cols : A.columns()-col);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
