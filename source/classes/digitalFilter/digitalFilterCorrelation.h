/***********************************************/
/**
* @file digitalFilterCorrelation.h
*
* @brief Simple correlation of time series.
*
* @author Matthias Ellmer
* @date 2017-11-07
*
*/
/***********************************************/

#ifndef __GROOPS_DIGITALFILTERCORRELATION__
#define __GROOPS_DIGITALFILTERCORRELATION__

// Latex documentation
#ifdef DOCSTRING_DigitalFilter
static const char *docstringDigitalFilterCorrelation = R"(
\subsection{Correlation}
Correlation ($\rho$) of \config{corr} is introduced into the time series:
\begin{equation}
  y_n = \rho\cdot y_{n-1} + \sqrt{1-\rho^2}x_n.
\end{equation}
)";
#endif

/***********************************************/

#include "classes/digitalFilter/digitalFilter.h"

/***** CLASS ***********************************/

/** @brief Correlation filter.
* @ingroup digitalFilterGroup
* @see DigitalFilter */
class DigitalFilterCorrelation : public DigitalFilterARMA
{
public:
  DigitalFilterCorrelation(Config &config);
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline DigitalFilterCorrelation::DigitalFilterCorrelation(Config &config)
{
  try
  {
    Double corr;
    readConfig(config, "correlation",       corr,              Config::MUSTSET, "",  "correlation");
    readConfig(config, "backwardDirection", backward,          Config::DEFAULT,  "0", "apply filter in backward direction");
    readConfig(config, "inFrequencyDomain", inFrequencyDomain, Config::DEFAULT,  "0", "apply filter in frequency domain");
    readConfig(config, "padType",           padType,           Config::MUSTSET,  "",  "");
    if(isCreateSchema(config)) return;

    bnStartIndex = 0;
    an = {1., -corr};
    bn = {std::sqrt(1-corr*corr)};
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
