/***********************************************/
/**
* @file matrixGeneratorSetType.h
*
* @brief Set type of a matrix.
*
* @author Torsten Mayer-Guerr
* @date 2017-09-01
*
*/
/***********************************************/

#ifndef __GROOPS_MATRIXGENERATORSETTYPE__
#define __GROOPS_MATRIXGENERATORSETTYPE__

// Latex documentation
#ifdef DOCSTRING_MatrixGenerator
static const char *docstringMatrixGeneratorSetType = R"(
\subsection{Set type}
Set type (matrix, matrixSymmetricUpper, matrixSymmetricLower, matrixTriangularUpper, matrixTriangularLower)
of a matrix. If the type is not matrix, the matrix must be quadratic. Symmetric matrices are filled symmetric
and for triangular matrix the other triangle is set to zero.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "matrixGenerator.h"

/***** CLASS ***********************************/

/** @brief Set type of a matrix.
* @ingroup matrixGeneratorGroup
* @see MatrixGenerator */
class MatrixGeneratorSetType : public MatrixGeneratorBase
{
  enum Type {GENERAL, SYMMETRICUPPER, SYMMETRICLOWER, TRIANGULARUPPER, TRIANGULARLOWER};
  Type type;
  MatrixGeneratorPtr matrix;

public:
  MatrixGeneratorSetType(Config &config);
  void compute(Matrix &A, UInt &startRow, UInt &startCol);
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline MatrixGeneratorSetType::MatrixGeneratorSetType(Config &config) : MatrixGeneratorBase(config)
{
  try
  {
    type = GENERAL;
    std::string choice;
    readConfig(config, "matrix", matrix, Config::MUSTSET, "", "");
    if(readConfigChoice(config, "type", choice, Config::MUSTSET, "", ""))
    {
      if(readConfigChoiceElement(config, "matrix",                choice, "")) type = GENERAL;
      if(readConfigChoiceElement(config, "matrixSymmetricUpper",  choice, "")) type = SYMMETRICUPPER;
      if(readConfigChoiceElement(config, "matrixSymmetricLower",  choice, "")) type = SYMMETRICLOWER;
      if(readConfigChoiceElement(config, "matrixTriangularUpper", choice, "")) type = TRIANGULARUPPER;
      if(readConfigChoiceElement(config, "matrixTriangularLower", choice, "")) type = TRIANGULARLOWER;
      endChoice(config);
    }
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void MatrixGeneratorSetType::compute(Matrix &A, UInt &/*startRow*/, UInt &/*startCol*/)
{
  try
  {
    A = matrix->compute();
    if((type != GENERAL) && (A.rows() != A.columns()))
      throw(Exception("Matrix must be quadratic"));
    if(type == GENERAL)         {A.setType(Matrix::GENERAL);}
    if(type == SYMMETRICUPPER)  {A.setType(Matrix::SYMMETRIC,  Matrix::UPPER); fillSymmetric(A);}
    if(type == SYMMETRICLOWER)  {A.setType(Matrix::SYMMETRIC,  Matrix::LOWER); fillSymmetric(A);}
    if(type == TRIANGULARUPPER) {A.setType(Matrix::TRIANGULAR, Matrix::UPPER); zeroUnusedTriangle(A);}
    if(type == TRIANGULARLOWER) {A.setType(Matrix::TRIANGULAR, Matrix::LOWER); zeroUnusedTriangle(A);}
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
