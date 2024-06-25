/***********************************************/
/**
* @file slrDesignMatrix.h
*
* @brief Management of sparse design matrix.
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#ifndef __GROOPS_SLRDESIGNMATRIX__
#define __GROOPS_SLRDESIGNMATRIX__

#include "parallel/matrixDistributed.h"
#include "slr/slrNormalEquationInfo.h"

/** @addtogroup slrGroup */
/// @{

/***** CLASS ***********************************/

/** @brief Management of sparse design matrix. */
class SlrDesignMatrix
{
  const SlrNormalEquationInfo   &normalEquationInfo;
  std::vector<UInt>              blockIndices;
  std::vector<UInt>              indexUsedBlock;
  std::vector<std::vector<UInt>> indexUsedParameter;
  std::vector<std::vector<UInt>> countUsedParameter;
  UInt                           row, rows;
  Matrix                         A;

public:
  Vector l;

  SlrDesignMatrix(const SlrNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef l=Vector());
 ~SlrDesignMatrix() {}

  void             init(const_MatrixSliceRef l);
  SlrDesignMatrix &selectRows(UInt row, UInt rows);
  MatrixSlice      column(const SlrParameterIndex &index);
  Matrix           mult(const_MatrixSliceRef x);
  Matrix           mult(const std::vector<Matrix> &x, UInt startBlock, UInt countBlock);
  void             transMult(const_MatrixSliceRef l, std::vector<Matrix> &x, UInt startBlock, UInt countBlock);
  void             accumulateNormals(MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount);
};

/// @}

/***********************************************/

#endif /* __GROOPS___ */
