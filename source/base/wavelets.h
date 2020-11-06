/***********************************************/
/**
* @file wavelets.h
*
* @brief Basic functions for wavelet decomposition
*
* @author Andreas Kvas
* @date 2017-06-28
*
*/
/***********************************************/

#ifndef __GROOPS_WAVELETS__
#define __GROOPS_WAVELETS__

#include "base/importStd.h"
#include "base/matrix.h"

/***** CLASS ***********************************/

/** @brief Wavelets.
* @ingroup base */
namespace Wavelets
{
  /** @brief FIR low-pass representation of wavelet. */
  Vector lowpass(const Vector &w);

  /** @brief FIR high-pass representation of wavelet. */
  Vector highpass(const Vector &w);

  /**
  * @brief Forward wavelet transform
  *
  * This routine computes the forward wavelet transform of a regularly sampled time series
  * using filter banks. The input time series is padded symmetrically around the edges to
  * reduce warmup effects. This implementation yields the same results as PyWavelets
  * (https://pywavelets.readthedocs.io/en/latest/).
  *
  * @param input data series
  * @param wl wavelet coefficients
  * @param maxLevel maximum decomposition level
  * @return standard vector containing each level of the decomposed time series  */
  std::vector<Matrix> waveletTransform(const Matrix &input, const Vector &wl, UInt maxLevel = MAX_UINT);

  /**
  * @brief Apply wavelet-halfband filters to time series and decimate.
  *
  * This routine applies the high- and low-pass filter representation of a wavelet
  * to a time series and decimates the result. The high-pass filtered input is stored in
  * detailCoefficients, the low-pass filtered input is stored in approxCoefficients.
  *
  * @param input data series
  * @param wl wavelet class object
  * @param[out] detailCoefficients high-pass filtered and decimated input
  * @param[out] approxCoefficients low-pass filtered and decimated input  */
  void halfbandFilter(const Matrix &input, const Vector &wl, Matrix &detailCoefficients, Matrix &approxCoefficients);

  /**
  * @brief Stationary wavelet transform
  *
  * This routine computes the stationary wavelet transform of a regularly sampled time series
  * using filter banks. As opposed to the regular wavelet transform, no decimation is performed, which means each level contains
  * the full number of samples. The input time series is padded symmetrically around the edges to
  * reduce warmup effects. This implementation yields the same results as PyWavelets (https://pywavelets.readthedocs.io/en/latest/).
  *
  * @param input data series
  * @param wl wavelet class object
  * @param maxLevel maximum decomposition level
  * @return standard vector containing each level of the decomposed time series
  */
  std::vector<Matrix> stationaryWaveletTransform(const Matrix &input, const Vector &wl, UInt maxLevel = MAX_UINT);

  /**
  * @brief Apply wavelet-halfband filters to time series.
  *
  * This routine applies the high- and low-pass filter representation of a wavelet
  * to a time series. As apposed to the standard half-band filter, no decimation is performed.
  * The high-pass filtered input is stored in detailCoefficients, the low-pass
  * filtered input is stored in approxCoefficients.
  *
  * @param input data series
  * @param wl wavelet class object
  * @param level level
  * @param[out] detailCoefficients high-pass filtered input
  * @param[out] approxCoefficients low-pass filtered input   */
  void stationaryHalfbandFilter(const Matrix &input, const Vector &wl, UInt level, Matrix &detailCoefficients, Matrix &approxCoefficients);

  /**
  * @brief Apply discrete wavelet transform to time series.
  *
  * This routine applies the discrete wavelet transform to time series
  * and stores the output as a vector with the sequence of approximation and details coefficients.
  *
  * @param signal data series
  * @param wl wavelet class object
  * @param J level
  * @param[out] coefficients wavelet coefficients
  * @param[out] flag padd length (if added) and level of decomposition
  * @param[out] length length of each set of details and approximation.*/
  void discreteWaveletTransform(std::vector<Double> &signal, const Vector &wl, UInt J, std::vector<Double> &coefficients, std::vector<UInt> &flag, std::vector<UInt> &length);

  /**
  * @brief Apply inverse discrete wavelet transform.
  *
  * This routine applies the inverse discrete wavelet transform to a set of coefficients.
  *
  * @param dwtop set of coefficients ${a_J, d_J,..., d_1}
  * @param wl wavelet class object
  * @param flag padd length (if added) and level of decomposition
  * @param recoveredSignal output in time domain
  * @param[out] length length of each set of details and approximation.*/
  void inverseDiscreteWaveletTransform(std::vector<Double> &dwtop, const Vector &wl, const std::vector<UInt> &flag, std::vector<Double> &recoveredSignal, std::vector<UInt> &length);

  /** @brief compute convolution of two vectors $a$ and $b$ and store the output in $c$.*/
  void convfft(std::vector<Double> &a, std::vector<Double> &b, std::vector<Double> &c);

  /** @brief compute fast fourier transform.*/
  void dfft(std::vector<std::complex<Double> > &data, int sign,UInt N);

  /** @brief apply a bit reverse algorithm to sort out the output frequencies of the FFT.*/
  void bitreverse(std::vector<std::complex<Double> > &signal);

  /** @brief compute decomposition and reconstruction filter coefficients from wavelet file.*/
  void filtCoef(const Vector &wl, std::vector<Double> &lpDecom, std::vector<Double> &hpDecom, std::vector<Double> &lpRecon, std::vector<Double> &hpRecon);

  /** @brief apply symmetric padding to a signal.*/
  void symmExtention(std::vector<Double> &sig, UInt paddLength);

  /** @brief apply down-sampling to a signal.*/
  void downSampling(std::vector<Double> &signal, UInt factor, std::vector<Double> &signalNew);

  /** @brief increase the sampling of a signal by inserting 0.*/
  void upSampling(std::vector<Double> &signal,   UInt factor, std::vector<Double> &signalNew);
}
#endif
