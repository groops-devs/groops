/***********************************************/
/**
* @file matrixGeneratorReshape.h
*
* @brief Matrix rechaped columnwise to new row and columns.
*
* @author Torsten Mayer-Guerr
* @date 2017-09-01
*
*/
/***********************************************/

#ifndef __GROOPS_MATRIXGENERATORRESHAPE__
#define __GROOPS_MATRIXGENERATORRESHAPE__

// Latex documentation
#ifdef DOCSTRING_MatrixGenerator
static const char *docstringMatrixGeneratorReshape = R"(
\subsection{Reshape}
Matrix reshaped columnwise to new row and columns.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "matrixGenerator.h"

/***** CLASS ***********************************/

/** @brief Matrix rechaped columnwise to new row and columns.
* @ingroup matrixGeneratorGroup
* @see MatrixGenerator */
class MatrixGeneratorReshape : public MatrixGeneratorBase
{
  ExpressionVariablePtr exprRows, exprCols;
  MatrixGeneratorPtr    matrix;

public:
  MatrixGeneratorReshape(Config &config);
  void compute(Matrix &A, UInt rowsBefore, UInt columnsBefore, UInt &startRow, UInt &startCol);
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline MatrixGeneratorReshape::MatrixGeneratorReshape(Config &config)
{
  try
  {
    readConfig(config, "matrix",  matrix,   Config::MUSTSET, "", "");
    readConfig(config, "rows",    exprRows, Config::MUSTSET, "0", "0: auto-determine rows, (variables: rowsBefore, columnsBefore)");
    readConfig(config, "columns", exprCols, Config::MUSTSET, "0", "0: auto-determine columns (variables: rowsBefore, columnsBefore)");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void MatrixGeneratorReshape::compute(Matrix &A, UInt rowsBefore, UInt columnsBefore, UInt &/*startRow*/, UInt &/*startCol*/)
{
  try
  {
    A = matrix->compute();

    VariableList varList;
    varList.setVariable("rowsBefore",    static_cast<Double>(rowsBefore));
    varList.setVariable("columnsBefore", static_cast<Double>(columnsBefore));
    varList.setVariable("rows",          static_cast<Double>(A.rows()));
    varList.setVariable("columns",       static_cast<Double>(A.columns()));
    A = reshape(A, static_cast<UInt>(exprRows->evaluate(varList)), static_cast<UInt>(exprCols->evaluate(varList)));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
