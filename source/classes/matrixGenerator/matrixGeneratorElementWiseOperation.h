/***********************************************/
/**
* @file matrixGeneratorElementWiseOperation.h
*
* @brief Elementwise operation between two matrices
*
* @author Andreas Kvas
* @date 2018-08-20
*
*/
/***********************************************/

#ifndef __GROOPS_MATRIXGENERATORELEMENTWISEOPERATION__
#define __GROOPS_MATRIXGENERATORELEMENTWISEOPERATION__

// Latex documentation
#ifdef DOCSTRING_MatrixGenerator
static const char *docstringMatrixGeneratorElementWiseOperation = R"(
\subsection{ElementWiseOperation}
Given two matrices $\mathbf{A}$ and $\mathbf{B}$ this class computes $c_{ij} = f(a_{ij}, b_{ij})$,
where $f$ is an expression (for example \verb|data0*data1|).
For each element of the matrix the variables \verb|data0|, \verb|data1|, \verb|row|, \verb|column|
are set and the expression \config{element} is evaluated.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "parser/dataVariables.h"
#include "matrixGenerator.h"

/***** CLASS ***********************************/

/** @brief Multiplication of matrices.
* @ingroup matrixGeneratorGroup
* @see MatrixGenerator */
class MatrixGeneratorElementWiseOperation : public MatrixGeneratorBase
{
  ExpressionVariablePtr expression;
  MatrixGeneratorPtr matrix1, matrix2;

public:
  MatrixGeneratorElementWiseOperation(Config &config);
  void compute(Matrix &A, UInt &startRow, UInt &startCol);
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline MatrixGeneratorElementWiseOperation::MatrixGeneratorElementWiseOperation(Config &config) : MatrixGeneratorBase(config)
{
  try
  {
    readConfig(config, "matrix1",    matrix1,    Config::MUSTSET, "",    "");
    readConfig(config, "matrix2",    matrix2,    Config::MUSTSET, "",    "");
    readConfig(config, "expression", expression, Config::DEFAULT, "data0*data1", "for each element of matrix (variables: data0, data1, row, column, rows, columns, rowsBefore, columnsBefore)");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void MatrixGeneratorElementWiseOperation::compute(Matrix &A, UInt &/*startRow*/, UInt &/*startCol*/)
{
  try
  {
    Matrix X = matrix1->compute();
    Matrix Y = matrix2->compute();
    if( (X.rows() != Y.rows()) || (X.columns() != Y.columns()) )
      throw(Exception("Matrix dimensions do not agree. ("+X.rows()%"%i x "s+X.columns()%"%i) vs. ("s+Y.rows()%"%i x "s+Y.columns()%"%i)."s));

    addVariable("rows",    static_cast<Double>(X.rows()),    varList);
    addVariable("columns", static_cast<Double>(X.columns()), varList);
    addVariable("row",     varList);
    addVariable("column",  varList);
    addVariable("data",    varList);
    expression->simplify(varList);

    A = Matrix(X.rows(), X.columns());
    for(UInt z=0; z<A.rows(); z++)
      for(UInt s=0; s<A.columns(); s++)
      {
        varList["row"]->setValue(static_cast<Double>(z));
        varList["column"]->setValue(static_cast<Double>(s));
        varList["data0"]->setValue(X(z,s));
        varList["data1"]->setValue(Y(z,s));
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
