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
  void compute(Matrix &A, UInt rowsBefore, UInt columnsBefore, UInt &startRow, UInt &startCol);
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline MatrixGeneratorExpression::MatrixGeneratorExpression(Config &config)
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

inline void MatrixGeneratorExpression::compute(Matrix &A, UInt rowsBefore, UInt columnsBefore, UInt &/*startRow*/, UInt &/*startCol*/)
{
  try
  {
    VariableList varList;
    varList.setVariable("rowsBefore",    static_cast<Double>(rowsBefore));
    varList.setVariable("columnsBefore", static_cast<Double>(columnsBefore));

    A = Matrix(static_cast<UInt>(exprRows->evaluate(varList)), static_cast<UInt>(exprCols->evaluate(varList)));

    varList.setVariable("rows",          static_cast<Double>(A.rows()));
    varList.setVariable("columns",       static_cast<Double>(A.columns()));
    varList.undefineVariable("row");
    varList.undefineVariable("column");
    expression->simplify(varList);

    for(UInt z=0; z<A.rows(); z++)
      for(UInt s=0; s<A.columns(); s++)
      {
        varList.setVariable("row",    static_cast<Double>(z));
        varList.setVariable("column", static_cast<Double>(s));
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
