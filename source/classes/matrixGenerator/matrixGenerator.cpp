/***********************************************/
/**
* @file matrixGenerator.cpp
*
* @brief Matrix calculation.
*
* @author Torsten Mayer-Guerr
* @date 2014-03-18
*
*/
/***********************************************/

#define DOCSTRING_MatrixGenerator

#include "base/import.h"
#include "config/configRegister.h"
#include "inputOutput/logging.h"
#include "classes/matrixGenerator/matrixGeneratorFile.h"
#include "classes/matrixGenerator/matrixGeneratorNormalsFile.h"
#include "classes/matrixGenerator/matrixGeneratorExpression.h"
#include "classes/matrixGenerator/matrixGeneratorElementManipulation.h"
#include "classes/matrixGenerator/matrixGeneratorElementWiseOperation.h"
#include "classes/matrixGenerator/matrixGeneratorAppend.h"
#include "classes/matrixGenerator/matrixGeneratorShift.h"
#include "classes/matrixGenerator/matrixGeneratorSlice.h"
#include "classes/matrixGenerator/matrixGeneratorReshape.h"
#include "classes/matrixGenerator/matrixGeneratorReorder.h"
#include "classes/matrixGenerator/matrixGeneratorSort.h"
#include "classes/matrixGenerator/matrixGeneratorTranspose.h"
#include "classes/matrixGenerator/matrixGeneratorMultiplication.h"
#include "classes/matrixGenerator/matrixGeneratorInverse.h"
#include "classes/matrixGenerator/matrixGeneratorCholesky.h"
#include "classes/matrixGenerator/matrixGeneratorRankKUpdate.h"
#include "classes/matrixGenerator/matrixGeneratorEigenValues.h"
#include "classes/matrixGenerator/matrixGeneratorDiagonal.h"
#include "classes/matrixGenerator/matrixGeneratorFromDiagonal.h"
#include "classes/matrixGenerator/matrixGeneratorSetType.h"
#include "classes/matrixGenerator/matrixGenerator.h"

/***********************************************/

GROOPS_REGISTER_CLASS(MatrixGenerator, "matrixGeneratorType",
                      MatrixGeneratorFile,
                      MatrixGeneratorNormalsFile,
                      MatrixGeneratorExpression,
                      MatrixGeneratorElementManipulation,
                      MatrixGeneratorElementWiseOperation,
                      MatrixGeneratorAppend,
                      MatrixGeneratorShift,
                      MatrixGeneratorSlice,
                      MatrixGeneratorReshape,
                      MatrixGeneratorReorder,
                      MatrixGeneratorSort,
                      MatrixGeneratorTranspose,
                      MatrixGeneratorMultiplication,
                      MatrixGeneratorInverse,
                      MatrixGeneratorCholesky,
                      MatrixGeneratorRankKUpdate,
                      MatrixGeneratorEigenValues,
                      MatrixGeneratorDiagonal,
                      MatrixGeneratorFromDiagonal,
                      MatrixGeneratorSetType)

GROOPS_READCONFIG_UNBOUNDED_CLASS(MatrixGenerator, "matrixGeneratorType")

/***********************************************/

MatrixGenerator::MatrixGenerator(Config &config, const std::string &name)
{
  try
  {
    std::string type;
    while(readConfigChoice(config, name, type, Config::OPTIONAL, "", "matrix calculation"))
    {
      if(readConfigChoiceElement(config, "file",                type, "from file"))
        matrix.push_back(new MatrixGeneratorFile(config));
      if(readConfigChoiceElement(config, "normalsFile",         type, "from normal equation file"))
        matrix.push_back(new MatrixGeneratorNormalsFile(config));
      if(readConfigChoiceElement(config, "expression",          type, "matrix filled by an expression"))
        matrix.push_back(new MatrixGeneratorExpression(config));
      if(readConfigChoiceElement(config, "elementManipulation", type, "elements of a matrix are manipulated by an expression"))
        matrix.push_back(new MatrixGeneratorElementManipulation(config));
      if(readConfigChoiceElement(config, "elementWiseOperation", type, "element wise operation of two matrices"))
        matrix.push_back(new MatrixGeneratorElementWiseOperation(config));
      if(readConfigChoiceElement(config, "append",              type, "append matrix to right or bottom"))
        matrix.push_back(new MatrixGeneratorAppend(config));
      if(readConfigChoiceElement(config, "shift",               type, "shift start row and/or start column"))
        matrix.push_back(new MatrixGeneratorShift(config));
      if(readConfigChoiceElement(config, "slice",               type, "slice of a matrix"))
        matrix.push_back(new MatrixGeneratorSlice(config));
      if(readConfigChoiceElement(config, "reshape",             type, "matrix reshaped columnwise to new row and columns"))
        matrix.push_back(new MatrixGeneratorReshape(config));
      if(readConfigChoiceElement(config, "reorder",             type, "reorder matrix with index vectors"))
        matrix.push_back(new MatrixGeneratorReorder(config));
      if(readConfigChoiceElement(config, "sort",                type, "sort matrix by column"))
        matrix.push_back(new MatrixGeneratorSort(config));
      if(readConfigChoiceElement(config, "transpose",           type, "transposed of a matrix"))
        matrix.push_back(new MatrixGeneratorTranspose(config));
      if(readConfigChoiceElement(config, "multiplication",      type, "multiplication of matrices"))
        matrix.push_back(new MatrixGeneratorMultiplication(config));
      if(readConfigChoiceElement(config, "inverse",             type, "inverse/pseudoinverse"))
        matrix.push_back(new MatrixGeneratorInverse(config));
      if(readConfigChoiceElement(config, "cholesky",            type, "cholesky decomposition"))
        matrix.push_back(new MatrixGeneratorCholesky(config));
      if(readConfigChoiceElement(config, "rankKUpdate",         type, "rank k update (A^T*A)"))
        matrix.push_back(new MatrixGeneratorRankKUpdate(config));
      if(readConfigChoiceElement(config, "eigenValues",         type, "eigen values"))
        matrix.push_back(new MatrixGeneratorEigenValues(config));
      if(readConfigChoiceElement(config, "diagonal",            type, "diagonal of matrix"))
        matrix.push_back(new MatrixGeneratorDiagonal(config));
      if(readConfigChoiceElement(config, "fromDiagonal",        type, "create matrix from diagonal vector"))
        matrix.push_back(new MatrixGeneratorFromDiagonal(config));
      if(readConfigChoiceElement(config, "setType",             type, "set type of a matrix"))
        matrix.push_back(new MatrixGeneratorSetType(config));
      endChoice(config);
      if(isCreateSchema(config))
        return;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

MatrixGenerator::~MatrixGenerator()
{
  for(UInt i=0; i<matrix.size(); i++)
    delete matrix.at(i);
}

/***********************************************/

Matrix MatrixGenerator::compute()
{
  try
  {
    std::vector<Matrix> A(matrix.size());
    std::vector<UInt>   startRow(matrix.size(), 0);
    std::vector<UInt>   startCol(matrix.size(), 0);
    UInt                rows = 0;
    UInt                cols = 0;
    for(UInt i=0; i<matrix.size(); i++)
    {
      matrix.at(i)->compute(A.at(i), rows, cols, startRow.at(i), startCol.at(i));
      rows = std::max(rows, A.at(i).rows()    + startRow.at(i));
      cols = std::max(cols, A.at(i).columns() + startCol.at(i));
    }

    // quick return?
    if((A.size() == 1) && (startRow.at(0) == 0) && (startCol.at(0) == 0))
      return A.at(0);

    // Accumulate all matrices
    Matrix B(rows, cols);
    Bool isTypeDefined = FALSE;
    for(UInt i=0; i<A.size(); i++)
      if(A.at(i).size())
      {
        // check type
        if(((A.at(i).getType() == Matrix::SYMMETRIC) || (A.at(i).getType() == Matrix::TRIANGULAR)) && (startRow.at(i) == startCol.at(i)))
        {
          if(!isTypeDefined)
            B.setType(A.at(i).getType(), A.at(i).isUpper() ? Matrix::UPPER : Matrix::LOWER);
          else if((A.at(i).getType() != B.getType()) || (A.at(i).isUpper() != B.isUpper()))
            A.at(i).setType(Matrix::GENERAL);
        }
        else if(A.at(i).getType() == Matrix::GENERAL)
          A.at(i).setType(Matrix::GENERAL);
        isTypeDefined = TRUE;

        // add
        axpy(1., A.at(i), B.slice(startRow.at(i), startCol.at(i), A.at(i).rows(), A.at(i).columns()));
      }

    return B;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
