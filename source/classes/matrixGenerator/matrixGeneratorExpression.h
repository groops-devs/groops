/***********************************************/
/**
* @file matrixGeneratorExpression.h
*
* @brief Matrix filled by an expression.
*
* @author Torsten Mayer-Guerr
* @date 2017-09-01
*
*/
/***********************************************/

#ifndef __GROOPS_MATRIXGENERATOREXPRESSION__
#define __GROOPS_MATRIXGENERATOREXPRESSION__

// Latex documentation
#ifdef DOCSTRING_MatrixGenerator
static const char *docstringMatrixGeneratorExpression = R"(
\subsection{Expression}
Matrix filled by an expression. For each element of the new matrix the variables
\verb|row| and \verb|column| are set and the expression \config{element} is evaluated.

Excample: The \config{element}=\verb|if(row==column,1,0)| generates an identity matrix.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "matrixGenerator.h"

/***** CLASS ***********************************/

/** @brief Matrix filled by an expression.
* @ingroup matrixGeneratorGroup
* @see MatrixGenerator */
class MatrixGeneratorExpression : public MatrixGeneratorBase
{
  ExpressionVariablePtr exprRows, exprCols;
  ExpressionVariablePtr expression;

public:
  MatrixGeneratorExpression(Config &config);
  void compute(Matrix &A, UInt &startRow, UInt &startCol);
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline MatrixGeneratorExpression::MatrixGeneratorExpression(Config &config) : MatrixGeneratorBase(config)
{
  try
  {
    readConfig(config, "rows",    exprRows,   Config::MUSTSET,  "", "(variables: rowsBefore, columnsBefore)");
    readConfig(config, "columns", exprCols,   Config::MUSTSET,  "", "(variables: rowsBefore, columnsBefore)");
    readConfig(config, "element", expression, Config::MUSTSET,  "if(row==column,1,0)",  "for each element of matrix (variables: row, column, rows, columns, rowsBefore, columnsBefore)");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void MatrixGeneratorExpression::compute(Matrix &A, UInt &/*startRow*/, UInt &/*startCol*/)
{
  try
  {
    A = Matrix(static_cast<UInt>(exprRows->evaluate(varList)), static_cast<UInt>(exprCols->evaluate(varList)));

    addVariable("rows",    static_cast<Double>(A.rows()),    varList);
    addVariable("columns", static_cast<Double>(A.columns()), varList);
    expression->simplify(varList);
    addVariable("row",    varList);
    addVariable("column", varList);

    for(UInt z=0; z<A.rows(); z++)
      for(UInt s=0; s<A.columns(); s++)
      {
        varList["row"]->setValue(static_cast<Double>(z));
        varList["column"]->setValue(static_cast<Double>(s));
        A(z,s) = expression->evaluate(varList);
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
