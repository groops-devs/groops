/***********************************************/
/**
* @file matrixGeneratorElementManipulation.h
*
* @brief The elements of a matrix are replaced an expression
*
* @author Torsten Mayer-Guerr
* @date 2017-09-01
*
*/
/***********************************************/

#ifndef __GROOPS_MATRIXGENERATORELEMENTMANIPULATION__
#define __GROOPS_MATRIXGENERATORELEMENTMANIPULATION__

// Latex documentation
#ifdef DOCSTRING_MatrixGenerator
static const char *docstringMatrixGeneratorElementManipulation = R"(
\subsection{Element manipulation}
The elements of a matrix are replaced an expression.
For each element of the matrix the variables \verb|data|, \verb|row|, \verb|column|
are set and the expression \config{element} is evaluated and replaces the element.
Additionally the standard data variables are available (assigned each row),
see~\reference{dataVariables}{general.parser:dataVariables}.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "parser/dataVariables.h"
#include "matrixGenerator.h"

/***** CLASS ***********************************/

/** @brief The elements of a matrix are replaced an expression
* @ingroup matrixGeneratorGroup
* @see MatrixGenerator */
class MatrixGeneratorElementManipulation : public MatrixGeneratorBase
{
  ExpressionVariablePtr expression;
  MatrixGeneratorPtr    matrix;

public:
  MatrixGeneratorElementManipulation(Config &config);
  void compute(Matrix &A, UInt rowsBefore, UInt columnsBefore, UInt &startRow, UInt &startCol);
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline MatrixGeneratorElementManipulation::MatrixGeneratorElementManipulation(Config &config)
{
  try
  {
    readConfig(config, "matrix",  matrix,     Config::MUSTSET, "",     "");
    readConfig(config, "element", expression, Config::MUSTSET, "data", "for each element of matrix (variables: data, row, column, rows, columns, rowsBefore, columnsBefore)");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void MatrixGeneratorElementManipulation::compute(Matrix &A, UInt rowsBefore, UInt columnsBefore, UInt &/*startRow*/, UInt &/*startCol*/)
{
  try
  {
    A = matrix->compute();

    VariableList varList;
    varList.setVariable("rowsBefore",    static_cast<Double>(rowsBefore));
    varList.setVariable("columnsBefore", static_cast<Double>(columnsBefore));
    varList.setVariable("rows",          static_cast<Double>(A.rows()));
    varList.setVariable("columns",       static_cast<Double>(A.columns()));

    addDataVariables(A, varList);
    varList.undefineVariable("row");
    varList.undefineVariable("column");
    varList.undefineVariable("data");
    expression->simplify(varList);

    for(UInt z=0; z<A.rows(); z++)
    {
      evaluateDataVariables(A, z, varList);
      for(UInt s=0; s<A.columns(); s++)
      {
        varList.setVariable("row",    static_cast<Double>(z));
        varList.setVariable("column", static_cast<Double>(s));
        varList.setVariable("data",   A(z,s));
        A(z,s) = expression->evaluate(varList);
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
