/***********************************************/
/**
* @file gnssLambda.h
*
* @brief LAMBDA - Least-squares AMBiguity Decorrelation Adjustment.
*
* @author Torsten Mayer-Guerr
* @date 2013-06-24
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSLAMBDA__
#define __GROOPS_GNSSLAMBDA__

/***********************************************/

#include <map>
#include "base/import.h"
#include "base/gnssType.h"

/***** CLASS ***********************************/

/** @brief LAMBDA - Least-squares AMBiguity Decorrelation Adjustment.
* @ingroup gnssGroup */
namespace GnssLambda
{
  /** @brief Sparse transformation matrix.
  * In both directions: direct and back.
  * Created by elementary transformations (swap and reduce) */
  class Transformation
  {
    std::vector<std::map<UInt, Double>> columnForward, columnBackward; // for each row a column with a factor

  public:
    Transformation(UInt dim=0);
   ~Transformation() {}

    void reduce(Double alpha, UInt i, UInt k); // row(k) -= alpha * row(i)
    void swap(UInt i, UInt k);

    Matrix transform(const_MatrixSliceRef x) const;
    Matrix transformBack(const_MatrixSliceRef x)  const;
    Matrix distributeBack(const_MatrixSliceRef x) const;
  };

  // ======================

  enum class IncompleteAction {STOP, SHRINKBLOCKSIZE, IGNORE, EXCEPTION};

  /** @brief Decorrelate ambiguities (Melbourne Wuebbena like linear combinations).
  * @param types list of phase observations.
  * @param wavelengthFactor 0.5 for old receivers using squaring technology.
  * @return Transformation matrix from decorrelated ambiguities to phase observations [cycles]->[m]. */
  Matrix phaseDecorrelation(const std::vector<GnssType> &types, Double wavelengthFactor);

  // LAMBDA method
  void   choleskyReversePivot(Matrix &N, Transformation &Z, UInt index0Z, Bool timing);
  Bool   choleskyReduce(UInt i, UInt k, MatrixSliceRef W, Transformation &transformation);
  Vector choleskyTransform(MatrixSliceRef W, Transformation &transformation, Bool timing);
  Bool   searchInteger(const_MatrixSliceRef xFloat, MatrixSliceRef W, const_MatrixSliceRef d, UInt maxSearchSteps, Vector &solution, Double &minNorm);
  Vector searchIntegerBlocked(const_MatrixSliceRef xFloat, MatrixSliceRef W, const_MatrixSliceRef d,
                              Double sigmaMaxResolve, UInt searchBlockSize, UInt maxSearchSteps, IncompleteAction incompleteAction, Bool timing,
                              Vector &isNotFixed, Double &sigma, Matrix &solutionSteps);
}

/***********************************************/

#endif /* __GROOPS___ */
