/***********************************************/
/**
* @file gnssDesignMatrix.h
*
* @brief Management of sparse design matrix.
*
* @author Torsten Mayer-Guerr
* @date 2018-04-01
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSDESIGNMATRIX__
#define __GROOPS_GNSSDESIGNMATRIX__

#include "parallel/matrixDistributed.h"
#include "gnss/gnssNormalEquationInfo.h"

/** @addtogroup gnssGroup */
/// @{

/***** CLASS ***********************************/

/** @brief Management of sparse design matrix. */
class GnssDesignMatrix
{
  const GnssNormalEquationInfo  &normalEquationInfo;
  std::vector<UInt>              blockIndices;
  std::vector<UInt>              indexUsedBlock;
  std::vector<std::vector<UInt>> indexUsedParameter;
  std::vector<std::vector<UInt>> countUsedParameter;
  UInt                           row, rows;
  Matrix                         A;

public:
  GnssDesignMatrix(const GnssNormalEquationInfo &normalEquationInfo, UInt rows=0);
 ~GnssDesignMatrix() {}

  void              init(UInt rows);
  GnssDesignMatrix &selectRows(UInt row, UInt rows);
  MatrixSlice       column(const GnssParameterIndex &index);
  MatrixSlice       column(UInt block, UInt col, UInt cols);
  Matrix            mult(const_MatrixSliceRef x);
  Matrix            mult(const std::vector<Matrix> &x, UInt startBlock, UInt countBlock);
  void              transMult(const_MatrixSliceRef l, std::vector<Matrix> &x, UInt startBlock, UInt countBlock);
  static void       accumulateNormals(const GnssDesignMatrix &A, const_MatrixSliceRef l, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount);
  static void       axpy(const std::vector<UInt> &rowInA, const std::vector<Double> &factors, const GnssDesignMatrix &B, GnssDesignMatrix &A);
};

/// @}

/***********************************************/

#endif /* __GROOPS___ */
