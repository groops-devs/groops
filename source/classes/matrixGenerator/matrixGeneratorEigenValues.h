/***********************************************/
/**
* @file matrixGeneratorEigenValues.h
*
* @brief EigenValues of a matrix.
*
* @author Torsten Mayer-Guerr
* @date 2017-09-01
*
*/
/***********************************************/

#ifndef __GROOPS_MATRIXGENERATOREIGENVALUES__
#define __GROOPS_MATRIXGENERATOREIGENVALUES__

// Latex documentation
#ifdef DOCSTRING_MatrixGenerator
static const char *docstringMatrixGeneratorEigenValues = R"(
\subsection{EigenValues}
Computes the eigenvalues of a square matrix and gives a vector of eigenvalues for symmetric matrices
or a matrix with 2 columns with real and imaginary parts in general case.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "matrixGenerator.h"

/***** CLASS ***********************************/

/** @brief EigenValues of a matrix.
* @ingroup matrixGeneratorGroup
* @see MatrixGenerator */
class MatrixGeneratorEigenValues : public MatrixGeneratorBase
{
  MatrixGeneratorPtr matrix;
  Bool eigenVectors;

public:
  MatrixGeneratorEigenValues(Config &config);
  void compute(Matrix &A, UInt rowsBefore, UInt columnsBefore, UInt &startRow, UInt &startCol);
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline MatrixGeneratorEigenValues::MatrixGeneratorEigenValues(Config &config)
{
  try
  {
    readConfig(config, "matrix",       matrix,       Config::MUSTSET,  "",  "");
    readConfig(config, "eigenVectors", eigenVectors, Config::DEFAULT,  "0", "return eigen vectors instead of eigen values");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void MatrixGeneratorEigenValues::compute(Matrix &A, UInt /*rowsBefore*/, UInt /*columnsBefore*/, UInt &/*startRow*/, UInt &/*startCol*/)
{
  try
  {
    A = matrix->compute().trans();
    Matrix eigen, VL, VR;
    if(A.getType() == Matrix::SYMMETRIC)
    {
      eigen = eigenValueDecomposition(A, eigenVectors);

      // reverse order
      if(eigenVectors)
        for(UInt k=0; k<A.columns()/2; k++)
        {
          Vector tmp(A.column(k));
          copy(A.column((A.columns()-1-k)), A.column(k));
          copy(tmp, A.column((A.columns()-1-k)));
        }
      else
      {
        A = Matrix(eigen.rows(), eigen.columns());
        for(UInt i=0; i<eigen.rows(); i++)
          copy(eigen.row(i), A.row((A.rows()-1-i)));
      }
    }
    else
    {
      eigen = eigenValueDecomposition(A, VL, VR, eigenVectors);
      if(!eigenVectors)
        A = eigen;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
