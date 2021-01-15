/***********************************************/
/**
* @file noiseGeneratorWhite.h
*
* @brief Simple white noise.
* @see NoiseGenerator
*
* @author Matthias Ellmer
* @date 2013-09-18
*
*/
/***********************************************/

#ifndef __GROOPS_NOISEGENERATORWHITE__
#define __GROOPS_NOISEGENERATORWHITE__

// Latex documentation
#ifdef DOCSTRING_NoiseGenerator
static const char *docstringNoiseGeneratorWhite = R"(
\subsection{White}
The noise is Gaussian with a standard deviation \config{sigma}.
The noise is computed via a pseudo random sequence with a start value given
by \config{initRandom}. The same value always yields the same sequence.
Be careful in \reference{parallel}{general.parallelization} mode
as all nodes generates the same pseudo random sequence.
If this value is set to zero a real random value is used as starting value.
)";
#endif

/***********************************************/

#include <random>
#include "parallel/parallel.h"
#include "classes/noiseGenerator/noiseGenerator.h"

/***** CLASS ***********************************/

/** @brief Simple white noise.
 * @ingroup noiseGeneratorGroup
 * @see NoiseGenerator */
class NoiseGeneratorWhite : public NoiseGeneratorBase
{
  std::mt19937_64 generator;
  Double          sigma;

public:
  NoiseGeneratorWhite(Config &config);
  Matrix noise(UInt samples, UInt series);
  Vector covarianceFunction(UInt length, Double sampling);
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline NoiseGeneratorWhite::NoiseGeneratorWhite(Config &config)
{
  try
  {
    UInt start;

    readConfig(config, "sigma",      sigma, Config::MUSTSET,  "1.0", "standard deviation");
    readConfig(config, "initRandom", start, Config::DEFAULT,  "0",   "start value for pseudo random sequence, 0: real random");
    if(isCreateSchema(config)) return;

    if(start)
      generator.seed(start);
    else
    {
      std::random_device randomDevice;
      generator.seed(randomDevice());
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Matrix NoiseGeneratorWhite::noise(UInt samples, UInt series)
{
  try
  {
    std::normal_distribution<Double> gaussian(0, sigma);
    Matrix e(samples, series);
    for(UInt i=0; i<e.rows(); i++)
      for(UInt k=0; k<e.columns(); k++)
        e(i,k) = gaussian(generator);
    return e;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Vector NoiseGeneratorWhite::covarianceFunction(UInt length, Double /*sampling*/)
{
  Vector covariance(length);
  covariance(0) = sigma*sigma;
  return covariance;
}

/***********************************************/

#endif /* __GROOPS_NOISEGENERATORWHITE__ */
