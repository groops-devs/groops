/***********************************************/
/**
* @file varianceComponentEstimation.h
*
* @brief Variance Component Estimation (VCE).
*
* @author Torsten Mayer-Guerr
* @date 2011-10-17
*
*/
/***********************************************/

#ifndef __GROOPS_VARIANCECOMPONENTESTIMATION__
#define __GROOPS_VARIANCECOMPONENTESTIMATION__

#include "base/import.h"
#include "inputOutput/fileName.h"

/***********************************************/

/** @brief Variance Component Estimation (VCE).
* @ingroup miscGroup */
namespace Vce
{
  /** @brief Random vector to estimate trace of a matrix
  * The number of @p columns increases the reliability. */
  Matrix monteCarlo(UInt rows, UInt columns);

  /** @brief Estimates the standardDeviation in case of otuliers.
  * The quadratic sum of residuals @p ePe and the @p redundancy are computed with downweigthed data. */
  Double standardDeviation(Double ePe, Double redundancy, Double huber, Double huberPower);

  /** @brief Robust least squares adjustment with multiple right hand sides.
  * Each @p countGroup (e.g. 3 for x,y,z) of observations are down weighted
  * if the estimated residuals are greater than @a huber times sigma0
  * @f[ s = \sqrt(e'Pe/r) > h\sigma_0. @f] */
  Matrix robustLeastSquares(const_MatrixSliceRef A, const_MatrixSliceRef l, UInt countGroup,
                            Double huber, Double huberPower, UInt maxIter, Vector &sigma);

  /** @brief Robust least squares adjustment with multiple right hand sides.
  * Each group starting at @p indexGroup and size `indexGroup.at(i+1)-indexGroup.at(i)`
  * (e.g. 3 for x,y,z) of observations are down weighted
  * if the estimated residuals are greater than @a huber times sigma0
  * @f[ s = \sqrt(e'Pe/r) > h\sigma_0. @f] */
  Matrix robustLeastSquares(const_MatrixSliceRef A, const_MatrixSliceRef l, const std::vector<UInt> &indexGroup,
                            Double huber, Double huberPower, UInt maxIter, Vector &sigma);

  /** @brief Compute the cos transformation matrix.
  * The matrix is normalized.
  * Applying this matrix two times results in the original input.
  * @code
  *   Matrix A(length, length);
  *   A(i,0) = 1/sqrt(2.*(length-1))
  *   A(i,k) = 2*cos(PI*i*k/(length-1))/sqrt(2.*(length-1)) // for k=1..n-2
  *   A(i,lenght-1) = (-1)**i/sqrt(2.*(length-1))
  * @endcode */
  Matrix cosTransform(UInt length);

  /** @brief Read covariance function from file or construct default function.
  * The covariance function must be saved as matrix.
  * The first column contains the time steps in seconds.
  * @param name file name (if empty a default cov function is contructed).
  * @param length the length of the cov functions (number of epochs).
  * @param columns number of expected covariance functions (e.g. 3 for orbits (along, cross, radial))
  * @param sampling in seconds.
  * @return columns number of covariance functions with length (plus sampling column). */
  Matrix readCovarianceFunction(const FileName &name, UInt length, UInt columns, Double sampling);

  /** @brief Compute redundancy matrix.
  * This function is called for each arc with decorrelated observation equation (@a We, @a WA, @a WB).
  * The covariance matrix @a W must be computed from @a CosTransform of the @a PSD.
  @code                                                                                *
  // covariance function for each axis.
  UInt countAxis = PSD.columns(); // e.g. for orbits (along, cross, radial).
  Matrix covFunc = CosTransform * PSD;

  // Compute the covariance matrix from the covariance functions
  Matrix Cov(index.size()*countAxis, Matrix::SYMMETRIC, Matrix::UPPER);
  for(UInt i=0; i*countAxis < Cov.rows(); i++)
    for(UInt k=i; k*countAxis < Cov.rows(); k++)
      for(UInt idxAxis=0; idxAxis < countAxis; idxAxis++)
        Cov(i*countAxis+idxAxis, k*countAxis+idxAxis) = cov(index.at(k)-index.at(i), idxAxis);
  @endcode
  The function returns the matrix of redundancies @a R and the weighted residuals @a WWe:
  @code
  WWe := Sigma^-1 e with Sigma = W^T W
  R   := Sigma^-1 - Sigma^-1 * (B A) * N^-1 * (B A)^T * Sigma^-1
        = Sigma^-1 - W^(-1)*Q1*Q1^T*W^T - W^(-1)*Q2*WA*N^(-1)*N^T*WA^T*Q2^T*W^T
  @endcode
  * @param W Cholesky decomposition of the covariance matrix.
  * @param We Decorrelated residuals.
  * @param WA Decorrelated design matrix.
  * @param WB Decorrelated design matrix.
  * @param[out] R Matrix of redundancies
  * @param[out] WWe Sigma^-1 e
  */
  void redundancy(const_MatrixSliceRef W, const_MatrixSliceRef We,
                  const_MatrixSliceRef WA, const_MatrixSliceRef WB,
                  Matrix &R, Vector &WWe);

  /** @brief Estimate frequency-wise variance factors for toeplitz covariance matrix.
  * This function is called for each arc with the matrix of redundancies @a R and the
  * weighted residuals @a WWe (see: computeRedundancy() ). This function computes the
  * squared sum of residuals and the redundancy for each frequency (row of @a ePe) and
  * component/axis (column).
  * @param R Matrix of redundancies.
  * @param WWe Weighted residuals Sigma^-1 e.
  * @param index Index of observations in the covariance function.
  * @param sigma Accuracy of the arc.
  * @param CosTransform To transform PSD to covariance function.
  * @param Psd Approximate PSD of the covariance function.
  * @param[in,out] ePe Updated squared sum of residuals.
  * @param[in,out] redundancy Updated for each axis (column) and frequency (row).
  * @param[in,out] ePeSum The sum of @a ePe other all frequencies is added.
  * @param[in,out] redundancySum The @a redundancy of ePe other all frequencies is added. */
  void psd(const_MatrixSliceRef R, const_MatrixSliceRef WWe,
           const std::vector<UInt> &index, Double sigma, const_MatrixSliceRef CosTransform, const_MatrixSliceRef Psd,
           MatrixSliceRef ePe, MatrixSliceRef redundancy, Double &ePeSum, Double &redundancySum);

  /** @brief Estimate variance factor for arbitrary covariance matrix.
  * This function is called for each arc with the matrix of redundancies @a R and the
  * weighted residuals @a WWe (see: computeRedundancy() ). This function computes the
  * squared sum of residuals and the redundancy for an arbitrary covariance matrix @a Cov.
  * @param R Matrix of redundancies.
  * @param WWe Weighted residuals Sigma^-1 e.
  * @param Cov Covariance matrix for which to determine variance factor.
  * @param[in,out] ePe Square sum of residuals due to @a C is added.
  * @param[in,out] redundancy Redundancy component due to @a C is added. */
  void matrix(const_MatrixSliceRef R, const_MatrixSliceRef WWe, const_MatrixSliceRef Cov,
              Double &ePe, Double &redundancy);

  /** @brief Updates the PSD from residuals and redundancies.
  * The PSD is updated by the factor ePe/redundancy for each axis (column) and frequency (row).
  * The maximum change factor or 1/factor of all amplitudes is returned in @a maxFactor. */
  void estimatePsd(MatrixSliceRef ePe, MatrixSliceRef redundancy, MatrixSliceRef Psd, Double &maxFactor, Bool jointZeroFrequency=TRUE);

  /** @brief Robust estimation of the mean sigma.
  * @a sigma contains the accuracy of each arc.
  * Only 50% of the median values are used to compute the mean value. */
  Double meanSigma(const Vector &sigma);
}

/***********************************************/

#endif /* __GROOPS__ */
