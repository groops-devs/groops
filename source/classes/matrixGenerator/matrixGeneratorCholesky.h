/***********************************************/
/**
* @file matrixGeneratorCholesky.h
*
* @brief Cholesky decomposition of a symmetric matrix.
*
* @author Torsten Mayer-Guerr
* @date 2017-09-01
*
*/
/***********************************************/

#ifndef __GROOPS_MATRIXGENERATORCHOLESKY__
#define __GROOPS_MATRIXGENERATORCHOLESKY__

// Latex documentation
#ifdef DOCSTRING_MatrixGenerator
static const char *docstringMatrixGeneratorCholesky = R"(
\subsection{Cholesky}
Upper triangular natrix of the cholesky decomposition of a symmetric matrix $\M A=\M W^T\M W$.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "matrixGenerator.h"

/***** CLASS ***********************************/

/** @brief Cholesky decomposition of a symmetric matrix.
* @ingroup matrixGeneratorGroup
* @see MatrixGenerator */
class MatrixGeneratorCholesky : public MatrixGeneratorBase
{
  MatrixGeneratorPtr matrix;

public:
  MatrixGeneratorCholesky(Config &config);
  void compute(Matrix &A, UInt rowsBefore, UInt columnsBefore, UInt &startRow, UInt &startCol);
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline MatrixGeneratorCholesky::MatrixGeneratorCholesky(Config &config)
{
  try
  {
    readConfig(config, "matrix", matrix, Config::MUSTSET, "", "");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void MatrixGeneratorCholesky::compute(Matrix &A, UInt /*rowsBefore*/, UInt /*columnsBefore*/, UInt &/*startRow*/, UInt &/*startCol*/)
{
  try
  {
    A = matrix->compute();
    cholesky(A);
    zeroUnusedTriangle(A);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
