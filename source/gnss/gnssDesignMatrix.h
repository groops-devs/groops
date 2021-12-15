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
  Vector l;

  GnssDesignMatrix(const GnssNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef l=Vector());
 ~GnssDesignMatrix() {}

  void              init(const_MatrixSliceRef l);
  GnssDesignMatrix &selectRows(UInt row, UInt rows);
  MatrixSlice       column(const GnssParameterIndex &index);
  Matrix            mult(const_MatrixSliceRef x);
  Matrix            mult(const std::vector<Matrix> &x, UInt startBlock, UInt countBlock);
  void              transMult(const_MatrixSliceRef l, std::vector<Matrix> &x, UInt startBlock, UInt countBlock);
  void              accumulateNormals(MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount);
};

/// @}

/***********************************************/

#endif /* __GROOPS___ */
