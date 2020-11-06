/***********************************************/
/**
* @file matrixGeneratorMultiplication.h
*
* @brief Multiplication of matrices.
*
* @author Torsten Mayer-Guerr
* @date 2017-09-01
*
*/
/***********************************************/

#ifndef __GROOPS_MATRIXGENERATORMULTIPLICATION__
#define __GROOPS_MATRIXGENERATORMULTIPLICATION__

// Latex documentation
#ifdef DOCSTRING_MatrixGenerator
static const char *docstringMatrixGeneratorMultiplication = R"(
\subsection{Multiplication}
Multiplication of matrices.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "matrixGenerator.h"

/***** CLASS ***********************************/

/** @brief Multiplication of matrices.
* @ingroup matrixGeneratorGroup
* @see MatrixGenerator */
class MatrixGeneratorMultiplication : public MatrixGeneratorBase
{
  Double factor;
  MatrixGeneratorPtr matrix1, matrix2;

public:
  MatrixGeneratorMultiplication(Config &config);
  void compute(Matrix &A, UInt &startRow, UInt &startCol);
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline MatrixGeneratorMultiplication::MatrixGeneratorMultiplication(Config &config) : MatrixGeneratorBase(config)
{
  try
  {
    readConfig(config, "matrix1", matrix1, Config::MUSTSET,  "",    "");
    readConfig(config, "matrix2", matrix2, Config::MUSTSET,  "",    "");
    readConfig(config, "factor",  factor,  Config::DEFAULT,  "1.0", "");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void MatrixGeneratorMultiplication::compute(Matrix &A, UInt &/*startRow*/, UInt &/*startCol*/)
{
  try
  {
    A = factor * matrix1->compute() * matrix2->compute();
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
