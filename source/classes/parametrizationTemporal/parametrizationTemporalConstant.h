/***********************************************/
/**
* @file parametrizationTemporalConstant.h
*
* @brief Constant in time.
* @see ParametrizationTemporal
*
* @author Torsten Mayer-Guerr
* @date 2011-12-28
*
*/
/***********************************************/

#ifndef __GROOPS_PARAMETRIZATIONTEMPORALCONSTANT__
#define __GROOPS_PARAMETRIZATIONTEMPORALCONSTANT__

// Latex documentation
#ifdef DOCSTRING_ParametrizationTemporal
static const char *docstringParametrizationTemporalConstant = R"(
\subsection{Constant}\label{parametrizationTemporalType:constant}
Represents a parameter being constant in time in each \config{interval}.
)";
#endif

/***********************************************/

#include "classes/parametrizationTemporal/parametrizationTemporal.h"

/***** CLASS ***********************************/

/** @brief Constant in time.
* @ingroup parametrizationTemporalGroup
* @see ParametrizationTemporal */
class ParametrizationTemporalConstant : public ParametrizationTemporalBase
{
  std::vector<Time> times;
  Bool              isInterval, includeLastTime;
  UInt              idxStart, idxEnd;

public:
  ParametrizationTemporalConstant(Config &config);

  Bool setInterval(const Time &timeStart, const Time &timeEnd, Bool estimatePerArc);
  UInt parameterCount() const {return (idxEnd-idxStart);}
  void factors(const Time &time, UInt startIndex, std::vector<UInt> &index, std::vector<Double> &factor) const;
  void parameterName(std::vector<ParameterName> &name) const;
};

/***********************************************/

#endif
