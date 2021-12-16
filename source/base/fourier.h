/***********************************************/
/**
* @file fourier.h
*
* @brief FFT-functions.
*
* @author Torsten Mayer-Guerr
* @author Andreas Kvas
* @date 2004-10-25
*/
/***********************************************/

#ifndef __GROOPS_FOURIER__
#define __GROOPS_FOURIER__

#include "base/importStd.h"
#include "base/matrix.h"

/***** CLASS ***********************************/

/**
* @brief FFT-functions.
* @ingroup base
*/
namespace Fourier
{
/***** Complex FFT *****************************/

  /** @brief Forward fourier transform of a real periodic sequence.
  * A data sequence with @f$j=0\ldots n-1@f$ elements is transformed with
  * @f[ y_k = \sum_{j=0}^{n-1} x_j e^{-2\pi jk/n}, @f]
  * where @f$k = 0 \ldots [n/2]@f$.
  * This transform is normalized since a call of fft followed by a call of synthesis gives the original result.
  * @param data data series
  * @return complex representation of the fourier transform */
  std::vector<std::complex<Double>> fft(const Vector &data);

  /** @brief Backward transform of complex fourier coefficients.
  *
  * If @a countEven=FALSE the size of the output data is @f$n=2m-1@f$ and computed by
  * @f[ x_j = \frac{1}{n} \sum_{k=-(m-1)}^{m-1} y_{|k|} e^{2\pi jk/n}, @f]
  * where @f$m@f$ is the number of complex fourier coefficients.
  *
  * If @a countEven=TRUE the size of the output data is @f$n=2m-2@f$ and computed by
  * @f[ x_j = \frac{1}{n} \sum_{k=-(m-2)}^{m-2} y_{|k|} e^{2\pi jk/n} + \frac{1}{n} (-1)^j re(y_{m-1}). @f]
  * The imaginary part of the last coefficient is ignored.
  *
  * @param F complex fourier coefficients
  * @param countEven size of output vector is even
  * @return data series */
  Vector synthesis(const std::vector<std::complex<Double>> &F, Bool countEven);

  /** @brief Frequency computation.
  * This function creates a frequency vector of half the length of an input
  * data vector. The output vector contains frequencies measured in cycles per time.
  * @param count length of input sequence
  * @param dt sampling interval
  * @return frequency vector */
  Vector frequencies(UInt count, Double dt);

  /** @brief Convert complex FFT coefficients to amplitude and phase.
  * @param F complex fourier coefficients
  * @param[out] amplitude amplitude spectrum
  * @param[out] phase phase spectrum */
  void complex2AmplitudePhase(const std::vector<std::complex<Double>> &F, Vector &amplitude, Vector &phase);

  /***** Covariance transformation ***************/

  /** @brief Computes the one sided PSD from a given covariance function.
  * A covariance function @f$c@f$with @f$j=0\ldots n-1@f$ elements is transformed to a one sided PSD
  * by the cosine transformation
  * @f[ psd_k = 2dt\left(c_0 + c_{n-1} (-1)^k + \sum_{j=1}^{n-2} 2 c_j \cos(\pi jk/(n-1))\right). @f]
  * The PSD of white noise with a variance of @f$\sigma^2@f$ is constant for all frequencies
  * with the power @f$A=2\delta t\sigma^2@f$.
  * @param cov covariance function
  * @param dt sampling interval
  * @return psd */
  Vector covariance2psd(const Vector &cov, Double dt);

  /** @brief Computes the covariance function from a one sided PSD.
  * A one sided PSD @f$p@f$ with @f$j=0\ldots n-1@f$ elements is transformed into a covariance function
  * by the cosine transformation
  * @f[ c_k = \frac{1}{4dt(n-1)}\left(p_0 + p_{n-1} (-1)^k + \sum_{j=1}^{n-2} 2 p_j \cos(\pi jk/(n-1))\right). @f]
  * The variance of a PSD with constant power @f$A@f$ for all frequencies (white noise)
  * is @f$\sigma^2=cov(0)=\frac{A}{2\delta t}@f$ .
  * @param psd power spectral density (PSD)
  * @param dt sampling interval
  * @return covariance function */
  Vector psd2covariance(const Vector &psd, Double dt);
}

/***********************************************/

#endif /* __GROOPS_FOURIER__ */
