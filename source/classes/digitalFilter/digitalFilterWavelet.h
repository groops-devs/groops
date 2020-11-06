/***********************************************/
/**
* @file digitalFilterWavelet.h
*
* @brief Filter representation of a wavelet.
*
* @author Andreas Kvas
* @date 2017-06-02
*
*/
/***********************************************/

#ifndef __GROOPS_DIGITALFILTERWAVELET__
#define __GROOPS_DIGITALFILTERWAVELET__

// Latex documentation
#ifdef DOCSTRING_DigitalFilter
static const char *docstringDigitalFilterWavelet = R"(
\subsection{Wavelet}
Filter representation of a wavelet.
)";
#endif

/***********************************************/

#include "base/wavelets.h"
#include "files/fileMatrix.h"
#include "classes/digitalFilter/digitalFilter.h"

/***** CLASS ***********************************/

/** @brief Filter representation of a wavelet.
* @ingroup digitalFilterGroup
* @see DigitalFilter */
class DigitalFilterWavelet : public DigitalFilterARMA
{
public:
  DigitalFilterWavelet(Config &config);
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline DigitalFilterWavelet::DigitalFilterWavelet(Config& config)
{
  try
  {
    FileName fileNameWavelet;
    Bool     islowpass = TRUE;
    UInt     level;

    readConfig(config, "inputfileWavelet", fileNameWavelet, Config::MUSTSET, "{groopsDataDir}/wavelets/", "wavelet coefficients");
    std::string choice;
    if(readConfigChoice(config, "type", choice, Config::MUSTSET, "", "filter type"))
    {
      if(readConfigChoiceElement(config, "lowpass",  choice)) islowpass = TRUE;
      if(readConfigChoiceElement(config, "highpass", choice)) islowpass = FALSE;
      endChoice(config);
    }
    readConfig(config, "level",             level,             Config::DEFAULT,  "0", "compute filter for specific decomposition level");
    readConfig(config, "backwardDirection", backward,          Config::DEFAULT,  "0", "apply filter in backward direction");
    readConfig(config, "inFrequencyDomain", inFrequencyDomain, Config::DEFAULT,  "0", "apply filter in frequency domain");
    readConfig(config, "padType",           padType,           Config::MUSTSET,  "",  "");
    if(isCreateSchema(config)) return;

    Vector wl;
    readFileMatrix(fileNameWavelet, wl);
    const Vector hn = (islowpass) ? Wavelets::lowpass(wl) : Wavelets::highpass(wl);
    const UInt upsamplingFactor = static_cast<UInt>(1)<<level;

    // set numerator and denominator polynomials
    an = Vector(1, 1.);
    bn = Vector(hn.rows() * upsamplingFactor);
    for(UInt k = 0; k < hn.rows(); k++)
      bn(k * upsamplingFactor) = hn(k);
    bnStartIndex = 0; // causal filter
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
