/***********************************************/
/**
* @file matrixGeneratorInverse.h
*
* @brief Inverse of a matrix.
*
* @author Torsten Mayer-Guerr
* @date 2017-09-01
*
*/
/***********************************************/

#ifndef __GROOPS_MATRIXGENERATORINVERSE__
#define __GROOPS_MATRIXGENERATORINVERSE__

// Latex documentation
#ifdef DOCSTRING_MatrixGenerator
static const char *docstringMatrixGeneratorInverse = R"(
\subsection{Inverse}
Inverse of a matrix $\M A^{-1}$.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "matrixGenerator.h"

/***** CLASS ***********************************/

/** @brief Inverse of a matrix.
* @ingroup matrixGeneratorGroup
* @see MatrixGenerator */
class MatrixGeneratorInverse : public MatrixGeneratorBase
{
  Bool               pseudo;
  MatrixGeneratorPtr matrix;

public:
  MatrixGeneratorInverse(Config &config);
  void compute(Matrix &A, UInt rowsBefore, UInt columnsBefore, UInt &startRow, UInt &startCol);
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline MatrixGeneratorInverse::MatrixGeneratorInverse(Config &config)
{
  try
  {
    readConfig(config, "matrix",        matrix, Config::MUSTSET, "",  "");
    readConfig(config, "pseudoInverse", pseudo, Config::DEFAULT, "0", "compute pseudo inverse instead of regular one");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void MatrixGeneratorInverse::compute(Matrix &A, UInt /*rowsBefore*/, UInt /*columnsBefore*/, UInt &/*startRow*/, UInt &/*startCol*/)
{
  try
  {
    A = matrix->compute();
    if(pseudo)
      A = pseudoInverse(A);
    else
      inverse(A);
    if(A.getType() == Matrix::SYMMETRIC)
      fillSymmetric(A);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
