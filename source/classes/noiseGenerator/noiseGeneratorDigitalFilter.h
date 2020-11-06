/***********************************************/
/**
* @file noiseGeneratorDigitalFilter.h
*
* @brief Generate noise using digital filter.
* @see NoiseGenerator
*
* @author Matthias Ellmer
* @date 2015-11-03
*
*/
/***********************************************/

#ifndef __GROOPS_NOISEGENERATORDIGITALFILTER__
#define __GROOPS_NOISEGENERATORDIGITALFILTER__

// Latex documentation
#ifdef DOCSTRING_NoiseGenerator
static const char *docstringNoiseGeneratorDigitalFilter = R"(
\subsection{Filter}
Generated noise \configClass{noise}{noiseGeneratorType} is
filtered by a \configClass{filter}{digitalFilterType}.
)";
#endif

/***********************************************/

#include "classes/digitalFilter/digitalFilter.h"
#include "classes/noiseGenerator/noiseGenerator.h"

/***** CLASS ***********************************/

/** @brief Generate noise from PSD.
 * @ingroup noiseGeneratorGroup
 * @see NoiseGenerator */
class NoiseGeneratorDigitalFilter : public NoiseGeneratorBase
{
  DigitalFilterPtr  filter;
  NoiseGeneratorPtr noisePtr;
  UInt              warmup, over;

public:
  NoiseGeneratorDigitalFilter(Config &config);
  Matrix noise(UInt samples, UInt series);
  Vector covarianceFunction(UInt length, Double sampling);
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline NoiseGeneratorDigitalFilter::NoiseGeneratorDigitalFilter(Config &config)
{
  try
  {
    readConfig(config, "filter",             filter,   Config::MUSTSET,  "",  "digital filter");
    readConfig(config, "noise",              noisePtr, Config::MUSTSET,  "",  "Basis noise");
    readConfig(config, "warmupEpochCount",   warmup,   Config::DEFAULT,  "0", "number of additional epochs at before start and after end");
    readConfig(config, "overSamplingFactor", over,     Config::DEFAULT,  "1", "noise with multiple higher sampling -> filter -> decimate");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Matrix NoiseGeneratorDigitalFilter::noise(UInt samples, UInt series)
{
  try
  {
    Matrix n = filter->filter(noisePtr->noise(over*samples+2*warmup, series));
    Matrix e(samples, series);
    for(UInt i=0; i<e.rows(); i++)
      copy(n.row(warmup+i*over), e.row(i));
    return e;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Vector NoiseGeneratorDigitalFilter::covarianceFunction(UInt length, Double sampling)
{
  try
  {
    auto frequencyResponse = filter->frequencyResponse(length*2);
    Vector psdBase = Fourier::covariance2psd(noisePtr->covarianceFunction(length, sampling).column(1), sampling);
    Vector psd(length);

    for(UInt i=0; i<length; i++)
      psd(i) = std::pow( std::abs(frequencyResponse.at(i)), 2) * psdBase(i);

    return Fourier::psd2covariance(psd, sampling);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif /* __GROOPS_NOISEGENERATORDIGITALFILTER__ */
