/***********************************************/
/**
* @file digitalFilterLag.h
*
* @brief Digital lag filter.
*
* @author Saniya Behzadpour
* @date 2018-03-13
*
*/
/***********************************************/

#ifndef __GROOPS_DIGITALFILTERLAG__
#define __GROOPS_DIGITALFILTERLAG__

// Latex documentation
#ifdef DOCSTRING_DigitalFilter
static const char *docstringDigitalFilterLag = R"(
\subsection{TimeLag}
Lag operator in digital filter representation.
)";
#endif

/***********************************************/

#include "classes/digitalFilter/digitalFilter.h"

/***** CLASS ***********************************/

/** @brief Digital lag filter.
* @ingroup digitalFilterGroup
* @see DigitalFilter */
class DigitalFilterLag : public DigitalFilterARMA
{
public:
  DigitalFilterLag(Config &config);
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline DigitalFilterLag::DigitalFilterLag(Config &config)
{
  try
  {
    Int lag;

    readConfig(config, "lag",               lag,               Config::MUSTSET,  "",  "lag epochs: 1 (lag); -1 (lead)");
    readConfig(config, "inFrequencyDomain", inFrequencyDomain, Config::DEFAULT,  "0", "apply filter in frequency domain");
    readConfig(config, "padType",           padType,           Config::MUSTSET,  "",  "");
    if(isCreateSchema(config)) return;

    an = Vector(1, 1.);
    bn = Vector(2*(std::abs(lag))+1);
    bn(std::abs(lag)+lag) = 1.0;
    bnStartIndex = bn.rows()/2;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
