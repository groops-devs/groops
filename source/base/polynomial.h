/***********************************************/
/**
* @file polynomial.h
*
* @brief Interpolation/derivation by polynomial.
*
* @author Torsten Mayer-Guerr
* @date 2017-05-27
*
*/
/***********************************************/

#ifndef __GROOPS_INTERPOLATION__
#define __GROOPS_INTERPOLATION__

#include "base/importStd.h"
#include "base/matrix.h"
#include "base/time.h"

/***** CLASS ***********************************/

/** @brief Interpolation/derivation by polynomial.
* @ingroup base */
class Polynomial
{
  Bool              throwException;
  UInt              degree;
  std::vector<Time> times;
  Double            sampling;
  Bool              isLeastSquares;
  Double            range, extrapolation;
  std::vector<Bool> isPrecomputed;
  Matrix            W;

public:
  /// Constructor
  Polynomial() {}

  /** @brief Constructor.
  * @see init(const std::vector<Time> &, UInt, Bool). */
  Polynomial(const std::vector<Time> &times, UInt degree, Bool throwException=TRUE) {init(times, degree, throwException);}

  /** @brief Initialize the interpolator.
  * @param times epochs of the input data.
  * @param degree of the polynomial.
  * @param throwException otherwise the non-interpolated values filled with NaN.
  * @param leastSquares use least squares fit polynomial instead of interpolation polynomial.
  * @param range interpolation is only allowed if all nodal points are within this interval or
  *              in case of least squares fit use all points within [t-range, t+range)
  *              (given in [seconds] or if negative given in the input median sampling).
  * @param extrapolation data will be extrapolated only if the nearest nodal point is less than this distance away
  *              (given in [seconds] or if negative given in the input median sampling).
  * @param margin to accelerate the interpolation check input times for equidistant sampling within margin [seconds].  */
  void init(const std::vector<Time> &times, UInt degree, Bool throwException, Bool leastSquares,
            Double range, Double extrapolation, Double margin=1e-5);

  /** @brief Simplified initialization.
  * Initialized with @a leastSquares=FALSE, @p range=-(degree+1.1), and @p Double extrapolation=-1.1
  * (Allow inter/extrapolation of one missing epoch).
  * @see init(const std::vector<Time> &, UInt, Bool, Bool, Double, Double, Double). */
  void init(const std::vector<Time> &times, UInt degree, Bool throwException=TRUE) {init(times, degree, throwException, FALSE/*leastSquares*/, -(degree+1.1), -1.1, 1e-5);}

  /** @brief Interpolate a matrix to new epochs.
  * @param timesNew output epochs of the returned matrix.
  * @param A input data of time series
  * @param rowsPerEpoch e.g. for @a A with positions (x,y,z) per epoch in separated rows.
  * @param derivative compute the kth derivative.
  * @param adjoint Apply the adjoint (transpose) operator (rows of A relates to timesNew, ouput rows relates to times).
  * @return Interpolated matrix with timesNew.size()*rowsPerEpoch rows. */
  Matrix interpolate(const std::vector<Time> &timesNew, const_MatrixSliceRef A, UInt rowsPerEpoch, UInt derivative, Bool adjoint=FALSE) const;

  /** @brief Interpolate a matrix to new epochs.
  * @see interpolate(const std::vector<Time> &, const_MatrixSliceRef, UInt, UInt) */
  Matrix interpolate(const std::vector<Time> &timesNew, const_MatrixSliceRef A, UInt rowsPerEpoch=1) const {return interpolate(timesNew, A, rowsPerEpoch, 0);}

  /** @brief Compute derivatives of a time series.
  * @see interpolate(const std::vector<Time> &, const_MatrixSliceRef, UInt, UInt) */
  Matrix derivative(const std::vector<Time> &timesNew, const_MatrixSliceRef A, UInt rowsPerEpoch=1) const {return interpolate(timesNew, A, rowsPerEpoch, 1);}

  /** @brief Compute 2nd derivatives of a time series.
  * @see interpolate(const std::vector<Time> &, const_MatrixSliceRef, UInt, UInt) */
  Matrix derivative2nd(const std::vector<Time> &timesNew, const_MatrixSliceRef A, UInt rowsPerEpoch=1) const {return interpolate(timesNew, A, rowsPerEpoch, 2);}
};

/***********************************************/

#endif /* __GROOPS_POLYNOMIAL__ */
