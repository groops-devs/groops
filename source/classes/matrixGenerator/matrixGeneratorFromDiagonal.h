/***********************************************/
/**
* @file matrixGeneratorFromDiagonal.h
*
* @brief Generate a matrix from a diagonal vector.
*
* @author Andreas Kvas
* @date 2018-09-01
*
*/
/***********************************************/

#ifndef __GROOPS_MATRIXGENERATORFROMDIAGONAL__
#define __GROOPS_MATRIXGENERATORFROMDIAGONAL__

// Latex documentation
#ifdef DOCSTRING_MatrixGenerator
static const char *docstringMatrixGeneratorFromDiagonal = R"(
\subsection{FromDiagonal}
Generate a matrix from a diagonal vector.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "matrixGenerator.h"

/***** CLASS ***********************************/

/** @brief  Generate a matrix from a diagonal vector.
* @ingroup matrixGeneratorGroup
* @see MatrixGenerator */
class MatrixGeneratorFromDiagonal : public MatrixGeneratorBase
{
  MatrixGeneratorPtr matrix;
  Int index;

public:
  MatrixGeneratorFromDiagonal(Config &config);
  void compute(Matrix &A, UInt &startRow, UInt &startCol);
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline MatrixGeneratorFromDiagonal::MatrixGeneratorFromDiagonal(Config &config) : MatrixGeneratorBase(config)
{
  try
  {
    readConfig(config, "matrix",   matrix, Config::MUSTSET, "",  "(nx1) or (1xn) diagonal vector");
    readConfig(config, "diagonal", index,  Config::DEFAULT, "0", "zero: main diagonal, positive: superdiagonal, negative: subdiagonal");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void MatrixGeneratorFromDiagonal::compute(Matrix &A, UInt &/*startRow*/, UInt &/*startCol*/)
{
  try
  {
    A = matrix->compute();
    if(std::min(A.rows(), A.columns()) != 1)
      throw(Exception("Input must be a (nx1) or (1xn) vector."));

    Vector d(flatten(A));
    const UInt outputSize = d.rows() + std::abs(index);

    A = Matrix(outputSize, outputSize);
    for(UInt k = 0; k<d.rows(); k++)
      A(k-std::min(0, index), k+std::max(0, index)) = d(k);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
