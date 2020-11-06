/***********************************************/
/**
* @file matrixGeneratorTranspose.h
*
* @brief Transpose of a matrix.
*
* @author Torsten Mayer-Guerr
* @date 2017-09-01
*
*/
/***********************************************/

#ifndef __GROOPS_MATRIXGENERATORTRANSPOSE__
#define __GROOPS_MATRIXGENERATORTRANSPOSE__

// Latex documentation
#ifdef DOCSTRING_MatrixGenerator
static const char *docstringMatrixGeneratorTranspose = R"(
\subsection{Transpose}
Transposed of a matrix $\M A^T$.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "matrixGenerator.h"

/***** CLASS ***********************************/

/** @brief Transposed of a matrix.
* @ingroup matrixGeneratorGroup
* @see MatrixGenerator */
class MatrixGeneratorTranspose : public MatrixGeneratorBase
{
  MatrixGeneratorPtr matrix;

public:
  MatrixGeneratorTranspose(Config &config);
  void compute(Matrix &A, UInt &startRow, UInt &startCol);
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline MatrixGeneratorTranspose::MatrixGeneratorTranspose(Config &config) : MatrixGeneratorBase(config)
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

inline void MatrixGeneratorTranspose::compute(Matrix &A, UInt &/*startRow*/, UInt &/*startCol*/)
{
  try
  {
    A = matrix->compute().trans();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
