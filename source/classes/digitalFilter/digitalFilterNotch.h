/***********************************************/
/**
* @file digitalFilterNotch.h
*
* @brief Notch filter.
*
* Implemented af Christian Siemes' dissertation, page 106.
*
* @author Andreas Kvas
* @date 2016-06-21
*
*/
/***********************************************/

#ifndef __GROOPS_DIGITALFILTERNOTCH__
#define __GROOPS_DIGITALFILTERNOTCH__

// Latex documentation
#ifdef DOCSTRING_DigitalFilter
static const char *docstringDigitalFilterNotch = R"(
\subsection{Notch}
Implemented after Christian Siemes' dissertation, page 106.

\fig{!hb}{0.6}{DigitalFilter_notch}{fig:DigitalFilterNotch}{Amplitude response of a notch filter of order three with default settings.}
)";
#endif

/***********************************************/

#include "classes/digitalFilter/digitalFilter.h"

/***** CLASS ***********************************/

/** @brief Butterworth filter.
* @ingroup digitalFilterGroup
* @see DigitalFilter */
class DigitalFilterNotch : public DigitalFilterARMA
{
public:
  DigitalFilterNotch(Config &config);
};

/***********************************************/

#endif
