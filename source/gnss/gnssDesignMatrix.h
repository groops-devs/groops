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

#include "gnss/gnss.h"

/** @addtogroup gnssGroup */
/// @{

/***** CLASS ***********************************/

/** @brief Management of sparse design matrix. */
class Gnss::DesignMatrix
{
  const NormalEquationInfo      &normalEquationInfo;
  std::vector<UInt>              blockIndices;
  std::vector<UInt>              indexUsedBlock;
  std::vector<std::vector<UInt>> indexUsedParameter;
  std::vector<std::vector<UInt>> countUsedParameter;
  UInt                           row, rows;
  Matrix                         A;

public:
  Vector l;

  DesignMatrix(const NormalEquationInfo &normalEquationInfo, const_MatrixSliceRef l=Vector());
 ~DesignMatrix() {}

  void          init(const_MatrixSliceRef l);
  DesignMatrix &selectRows(UInt row_, UInt rows_);
  MatrixSlice   column(const Gnss::ParameterIndex &index);
  Matrix        mult(const_MatrixSliceRef x);
  Matrix        mult(const std::vector<Matrix> &x, UInt startBlock, UInt countBlock);
  void          transMult(const_MatrixSliceRef l, std::vector<Matrix> &x, UInt startBlock, UInt countBlock);
  void          accumulateNormals(MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount);
};

/// @}

/***********************************************/

#endif /* __GROOPS___ */
