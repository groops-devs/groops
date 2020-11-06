/***********************************************/
/**
* @file digitalFilterMovingAverage.h
*
* @brief Moving average (boxcar) filter.
*
* @author Andreas Kvas
* @date 2016-06-21
*
*/
/***********************************************/

#ifndef __GROOPS_DIGITALFILTERMOVINGAVERAGE__
#define __GROOPS_DIGITALFILTERMOVINGAVERAGE__

// Latex documentation
#ifdef DOCSTRING_DigitalFilter
static const char *docstringDigitalFilterMovingAverage = R"(
\subsection{MovingAverage}
Moving average (boxcar) filter. For odd lengths, this filter is symmetric and has therefore no phase shift. For even lengths, a phase shift of half a cycle is introduced.

\[
  y_n = \sum_{k=-\lfloor\frac{P}{2}\rfloor}^{\lfloor\frac{P}{2}\rfloor} \frac{1}{P}x_{n-k}
\]

)";
#endif

/***********************************************/

#include "classes/digitalFilter/digitalFilter.h"

/***** CLASS ***********************************/

/** @brief Moving average (boxcar) filter.
* @ingroup digitalFilterGroup
* @see DigitalFilter */
class DigitalFilterMovingAverage : public DigitalFilterARMA
{
public:
  DigitalFilterMovingAverage(Config &config);
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline DigitalFilterMovingAverage::DigitalFilterMovingAverage(Config &config)
{
  try
  {
    UInt length;

    readConfig(config, "length",            length,            Config::MUSTSET, "", "number of epochs in averaging operator");
    readConfig(config, "inFrequencyDomain", inFrequencyDomain, Config::DEFAULT,  "0", "apply filter in frequency domain");
    readConfig(config, "padType",           padType,           Config::MUSTSET, "", "");
    if(isCreateSchema(config)) return;

    an = Vector(1, 1.);
    bn = Vector(length);
    for(UInt k = 0; k<bn.rows(); k++)
      bn(k) = 1./length;
    bnStartIndex = bn.rows()/2;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
