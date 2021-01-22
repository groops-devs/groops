/***********************************************/
/**
* @file normalsShortTimeStaticLongTime.h
*
* @brief Normal equations with short and long time gravity variations.
*
* @author Torsten Mayer-Guerr
* @date 2020-11-29
*
*/
/***********************************************/

#ifndef __GROOPS_NORMALSSHORTTIMESTATICLONGTIME__
#define __GROOPS_NORMALSSHORTTIMESTATICLONGTIME__

#include "base/import.h"
#include "parallel/matrixDistributed.h"
#include "classes/observation/observation.h"

/***** TYPES ***********************************/

class ParameterSelector;
class ParametrizationTemporal;
typedef std::shared_ptr<ParameterSelector>       ParameterSelectorPtr;
typedef std::shared_ptr<ParametrizationTemporal> ParametrizationTemporalPtr;

/***** CLASS ***********************************/

/** @brief Normal equations with short and long time gravity variations.
* @ingroup miscGroup */
class NormalsShortTimeStaticLongTime : public MatrixDistributed
{
public:
  Matrix n;         // right hand sides
  Vector lPl;       // =l'Pl, weighted norm of the observations
  UInt   obsCount;  // number of observations

  std::vector<ParameterName>       parameterNames;
  std::vector<std::vector<UInt>>   indexA; // for each interval
  std::vector<std::vector<UInt>>   indexN; // for each interval
  UInt                             blockIndexStatic;
  std::vector<UInt>                blockIndexShortTime;
  UInt                             blockCountTemporal;
  std::vector<UInt>                blockIndexTemporal; // static + for each temporal
  std::vector<MatrixDistributed>   normalsTemporal;
  std::vector<Matrix>              nTemporal;
  std::vector<std::vector<Double>> factorTemporal;

  void init(ObservationPtr observation, const std::vector<Time> &timesInterval,
            UInt defaultBlockSize, Parallel::CommunicatorPtr comm, Bool sortStateBeforeGravityParameter,
            UInt countShortTimeParameters, ParameterSelectorPtr parameterShortTime,
            ParametrizationTemporalPtr temporal=nullptr, ParameterSelectorPtr parameterTemporal=nullptr);

  /// Allocate memory of matrix blocks.
  void setBlocks(const std::vector<UInt> &arcsInterval);

  /// Fill all matrix blocks with zero.
  void setNull();

  /** @brief Add observation equations to the normal system.
  * The design matrix @p A must contain only parameters valid in interval @p idInterval
  * (Computed with @a setInterval()).
  * The input matrices @p l, @p A, and @p B might be destroyed. */
  void accumulate(UInt idInterval, Matrix &l, Matrix &A, Matrix &B);

  /// Reduce the system of normal equations.
  void reduceSum(Bool timing=TRUE);

  void addShortTimeNormals(Double sigma2, const std::vector<std::vector<std::vector<Matrix>>> &normalsShortTime);

  void regularizeUnusedParameters(UInt countBlock);

  Double solve(Matrix &x, Matrix &Wz);

  /** @brief Compute the standard devitiations of the parameter vector.
  * (The sqrt of the diagonals of inverse normal matrix). The normals must be given as cholesky decomposition (with @a solve).
  * On oputput the normal matrix contains the inverse of cholesky matrix. */
  Vector parameterStandardDeviation();

  Double estimateShortTimeNormalsVariance(Double sigma2, const std::vector<std::vector<std::vector<Matrix>>> &normalsShortTime,
                                          const_MatrixSliceRef x, const_MatrixSliceRef Wz) const;
  /** @brief Computes Ax = factor*A*x.
  * The design matrix @p A contains only parameters valid in interval @p idInterval. */
  void designMatMult(UInt idInterval, Double factor, const_MatrixSliceRef A, const_MatrixSliceRef x, MatrixSliceRef Ax);
};

/***********************************************/

#endif
