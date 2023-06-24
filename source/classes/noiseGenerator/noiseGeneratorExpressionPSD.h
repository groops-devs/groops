/***********************************************/
/**
* @file noiseGeneratorExpressionPSD.h
*
* @brief Generate noise defined by one sided PSD.
* @see NoiseGenerator
*
* @author Torsten Mayer-Guerr
* @date 2017-09-07
*
*/
/***********************************************/

#ifndef __GROOPS_NOISEGENERATOREXPRESSIONPSD__
#define __GROOPS_NOISEGENERATOREXPRESSIONPSD__

// Latex documentation
#ifdef DOCSTRING_NoiseGenerator
static const char *docstringNoiseGeneratorExpressionPSD = R"(
\subsection{ExpressionPSD}
This generator creates noise defined by a one sided PSD.
The \config{psd} is an expression controlled by the variable 'freq'.
To determine the frequency \config{sampling} must be given.
)";
#endif

/***********************************************/

#include "base/fourier.h"
#include "inputOutput/logging.h"
#include "classes/noiseGenerator/noiseGenerator.h"

/***** CLASS ***********************************/

/** @brief Generate noise defined by one sided PSD.
 * @ingroup noiseGeneratorGroup
 * @see NoiseGenerator */
class NoiseGeneratorExpressionPSD : public NoiseGeneratorBase
{
  NoiseGeneratorPtr     noisePtr;
  ExpressionVariablePtr expression;
  Double                sampling;

public:
  NoiseGeneratorExpressionPSD(Config &config);
  Matrix noise(UInt samples, UInt series);
  Vector covarianceFunction(UInt length, Double sampling);
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline NoiseGeneratorExpressionPSD::NoiseGeneratorExpressionPSD(Config &config)
{
  try
  {
    readConfig(config, "noise",    noisePtr,   Config::MUSTSET, "",  "Basis noise");
    readConfig(config, "psd",      expression, Config::MUSTSET, "1", "one sided PSD (variable: freq [Hz]) [unit^2/Hz]");
    readConfig(config, "sampling", sampling,   Config::MUSTSET, "1", "to determine frequency [seconds]");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Matrix NoiseGeneratorExpressionPSD::noise(UInt samples, UInt series)
{
  try
  {
    // Fill the sequence wk with white Noise, second half needs to be zero for fourier padding
    Matrix wk(2*samples, series);
    copy(noisePtr->noise(samples, series), wk.row(0, samples));

    Vector freq = Fourier::frequencies(wk.rows(), sampling);
    Vector PSD(freq.rows());
    VariableList varList;
    for(UInt i=0; i<freq.rows(); i++)
    {
      varList.setVariable("freq", freq(i));
      PSD(i) = expression->evaluate(varList);
      if(!std::isfinite(PSD(i)))
        logWarning << "Warning: PSD at frequency "+freq(i) % "%f [Hz] is "s + PSD(i) % "%f"s << Log::endl;
    }

    // For each column, perform discrete Fourier transform, multiply the results and synthesize noise
    Matrix noise(samples, series);
    for(UInt k=0; k<series; k++)
    {
      auto Wk = Fourier::fft(Vector(wk.column(k))); // Complex Fourier transform
      for(UInt i=0; i<PSD.rows(); i++)
        Wk.at(i) *= std::sqrt(PSD(i));
      copy(Fourier::synthesis(Wk, TRUE/*even*/).slice(0, samples), noise.column(k));
    }

    return noise;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Vector NoiseGeneratorExpressionPSD::covarianceFunction(UInt length, Double sampling)
{

  Vector freq = Fourier::frequencies(2*length, sampling);
  Vector psd  = Fourier::covariance2psd(noisePtr->covarianceFunction(length, sampling).column(1), sampling);

  VariableList varList;
  for(UInt i=0; i<length; i++)
  {
    varList.setVariable("freq", freq(i));
    psd(i) *= expression->evaluate(varList);
  }

  return Fourier::psd2covariance(psd, sampling);
}

/***********************************************/
#endif /* __GROOPS_NOISEGENERATOR__ */
