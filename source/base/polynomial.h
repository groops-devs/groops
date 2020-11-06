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
  Matrix W;

public:
  /// Constructor
  Polynomial() {}

  /// Constructor
  explicit Polynomial(UInt degree) {init(degree);}

  /// Set the degree of the interpolation polynomial. */
  void init(UInt degree);

  /** @brief Interpolate a matrix to new epochs.
  * @param timesNew output epochs of the returned matrix.
  * @param times epoch of the rows of @a A (must be constant sampling)
  * @param A input data of time series
  * @param rowsPerEpoch e.g. for @a A with positions (x,y,z) per epoch in separated rows.
  * @return Interpolated matrix with timesNew.size()*rowsPerEpoch rows. */
  Matrix interpolate(const std::vector<Time> &timesNew, const std::vector<Time> &times, const_MatrixSliceRef A, UInt rowsPerEpoch=1) const;

  /** @brief Compute derivatives of a time series.
  * @param sampling in seconds
  * @param A input data of time series with constant sampling
  * @param rowsPerEpoch e.g. for @a A with positions (x,y,z) per epoch in separated rows.
  * @return Matrix with derivatives at same epochs as @a A. */
  Matrix derivative(Double sampling, const_MatrixSliceRef A, UInt rowsPerEpoch=1) const;

  /** @brief Compute 2nd derivatives of a time series.
  * @param sampling in seconds
  * @param A input data of time series with constant sampling
  * @param rowsPerEpoch e.g. for @a A with positions (x,y,z) per epoch in separated rows.
  * @return Matrix with 2nd derivatives at same epochs as @a A. */
  Matrix derivative2nd(Double sampling, const_MatrixSliceRef A, UInt rowsPerEpoch=1) const;

  /** @brief Integration of a time series.
  * @param sampling in seconds
  * @param A input data of time series with constant sampling
  * @param rowsPerEpoch e.g. for @a A with positions (x,y,z) per epoch in separated rows.
  * @return Integrals from start to epoch at each row (first row is zero). */
  Matrix integration(Double sampling, const_MatrixSliceRef A, UInt rowsPerEpoch=1) const;
};

/***********************************************/

#endif /* __GROOPS_POLYNOMIAL__ */
