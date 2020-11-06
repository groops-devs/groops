/***********************************************/
/**
* @file noiseGeneratorPowerLaw.h
*
* @brief Generate noise following a power law relationship like 1/f^alpha.
* @see NoiseGenerator
*
* @see Kasdin (1995): Discrete simulation of colored noise and stochastic processes and 1/f^alpha power law noise generation, Proceedings of the IEEE, Volume 83, Number 5, 1995, pages 802-827.
* @see Kasdin (1992): Discrete simulation of power law noise [for oscillator stability evaluation]
*
* @author Matthias Ellmer
* @date 2014-01-29
*
*/
/***********************************************/

#ifndef __GROOPS_NOISEGENERATORPOWERLAW__
#define __GROOPS_NOISEGENERATORPOWERLAW__

// Latex documentation
#ifdef DOCSTRING_NoiseGenerator
static const char *docstringNoiseGeneratorPowerLaw = R"(
\subsection{PowerLaw}
This generator creates noise that conforms to a power law relationship, where the power
of the noise at a frequency is proportional to $1/f^\alpha$, with a typically between -2 and 2.
)";
#endif

/***********************************************/

#include "base/fourier.h"
#include "classes/noiseGenerator/noiseGenerator.h"

/***** CLASS ***********************************/

/** @brief Generate noise following power law.
 * @ingroup noiseGeneratorGroup
 * @see NoiseGenerator */
class NoiseGeneratorPowerLaw : public NoiseGeneratorBase
{
  NoiseGeneratorPtr noisePtr;
  Double alpha;

public:
  NoiseGeneratorPowerLaw(Config &config);
  Matrix noise(UInt samples, UInt series);
  Vector covarianceFunction(UInt length, Double sampling);
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline NoiseGeneratorPowerLaw::NoiseGeneratorPowerLaw(Config &config)
{
  try
  {
    readConfig(config, "noise", noisePtr, Config::MUSTSET, "",  "Basis noise");
    readConfig(config, "alpha", alpha,    Config::MUSTSET, "0", "Exponent of the power law relationship 1/f^alpha");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

// This is a reimplemantation of the code given in [Kasdin,1995], Appendix 2
inline Matrix NoiseGeneratorPowerLaw::noise(UInt samples, UInt series)
{
  try
  {
    // Fill the sequence wk with white Noise, second half padded with zeros
    Matrix wk(2*samples, series);
    copy(noisePtr->noise(samples, series), wk.row(0, samples));

    // Generate coefficients hk in first half of vector, second half padded with zeros
    Vector hk(2*samples);
    hk(0) = 1.0;
    for(UInt i=1; i<samples; i++)
      hk(i) = hk(i-1) * (0.5*alpha+(i-1))/i;
    auto Hk = Fourier::fft(hk);

    // For each column, perform discrete Fourier transform, multiply the results and synthesize noise
    Matrix noise(samples, series);
    for(UInt i=0; i<series; i++)
    {
      auto Wk = Fourier::fft(wk.column(i));
      for(UInt j=0; j<Wk.size(); j++)
        Wk.at(j) *= Hk.at(j); // complex multiplication
      copy(Fourier::synthesis(Wk, TRUE/*even*/).slice(0, samples), noise.column(i));
    }

    return noise;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector NoiseGeneratorPowerLaw::covarianceFunction(UInt length, Double sampling)
{
  try
  {
    Vector freq = Fourier::frequencies(2*length, sampling);
    Vector psd  = Fourier::covariance2psd(noisePtr->covarianceFunction(length, sampling).column(1), sampling);

    // Skip zero frequency
    // eq. 98 in [Kasdin,1995]
    // I think it should be psd(i) = std::pow( freq(i) * (sigma*sigma), -alpha). Why it is not, I don't quite understand.
    for(UInt i=1; i<length; i++)
      psd(i) *= std::pow(2 * std::sin(PI * freq(i) * sampling), -alpha);

    return Fourier::psd2covariance(psd, sampling);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif /* __GROOPS_NOISEGENERATORPOWERLAW__ */
