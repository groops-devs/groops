/***********************************************/
/**
* @file digitalFilterButterworth.h
*
* @brief Butterworth filter.
*
* @author Matthias Ellmer
* @author Andreas Kvas
* @date 2016-06-21
*
*/
/***********************************************/

#ifndef __GROOPS_DIGITALFILTERBUTTERWORTH__
#define __GROOPS_DIGITALFILTERBUTTERWORTH__

// Latex documentation
#ifdef DOCSTRING_DigitalFilter
static const char *docstringDigitalFilterButterworth = R"(
\subsection{Butterworth}
Digital implementation of the Butterworth filter. The design of the filter is done by modifying the analog (continuous time) transfer function, which is
then transformed into the digital domain by using the bilinear transform. The filter coefficients are then determined by a least squares adjustment in time domain.

The \config{filterType} can be \config{lowpass}, \config{highpass}, where one cutoff frequency has to be specified, and \config{bandpass} and \config{bandstop} where to cutoff frequencies have to be specified.
Cutoff frequencies must be given as normalized frequency $w_n = f/f_{\text{nyq}}$. For a cutoff frequency of 30~mHz for a time series sampled with 5~seconds gives a normalized frequency of $0.03/0.1 = 0.3$.
)";
#endif

/***********************************************/

#include "classes/digitalFilter/digitalFilter.h"

/***** CLASS ***********************************/

/** @brief Butterworth filter.
* @ingroup digitalFilterGroup
* @see DigitalFilter */
class DigitalFilterButterworth : public DigitalFilterARMA
{
  static std::vector< std::complex<Double> > transferFunction(const std::vector<std::complex<Double> > &z, std::function<std::vector<std::complex<Double>>(const std::vector<std::complex<Double>>&)>& mapping, UInt order, Double g0=1.0);
  static std::vector< std::complex<Double> > butterworthPolynomial(const std::vector< std::complex<Double> > &s, UInt order);
  static std::vector< std::complex<Double> > lowpassFrequencyScaling (const std::vector< std::complex<Double> > &z, Double Wn);
  static std::vector< std::complex<Double> > highpassFrequencyScaling(const std::vector< std::complex<Double> > &z, Double Wn);
  static std::vector< std::complex<Double> > bandpassFrequencyScaling(const std::vector< std::complex<Double> > &z, Double Wn, Double Vn);
  static std::vector< std::complex<Double> > bandstopFrequencyScaling(const std::vector< std::complex<Double> > &z, Double Wn, Double Vn);

public:
  DigitalFilterButterworth(Config &config);
};

/***********************************************/

#endif
