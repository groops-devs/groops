/***********************************************/
/**
* @file parametrizationTemporalPolynomial.h
*
* @brief Legendre Polynomials in time intervals.
* @see ParametrizationTemporal
*
* @author Torsten Mayer-Guerr
* @date 2015-05-31
*
*/
/***********************************************/

#ifndef __GROOPS_PARAMETRIZATIONTEMPORALPOLYNOMIAL__
#define __GROOPS_PARAMETRIZATIONTEMPORALPOLYNOMIAL__

// Latex documentation
#ifdef DOCSTRING_ParametrizationTemporal
static const char *docstringParametrizationTemporalPolynomial = R"(
\subsection{Polynomial}
A time variable function is represented by Legendre polynomials in each \config{interval}.
The time is normed to $[-1,1)$ in each interval.

The total parameter count is $(N+1)M$,
where $N$ is the polynmial degree and $M$ the number of intervals.
)";
#endif

/***********************************************/

#include "classes/parametrizationTemporal/parametrizationTemporal.h"

/***** CLASS ***********************************/

/** @brief Splines in time domain.
* @ingroup parametrizationTemporalGroup
* @see ParametrizationTemporal */
class ParametrizationTemporalPolynomial : public ParametrizationTemporalBase
{
  std::vector<Time> times;
  Bool              isInterval;
  UInt              idxStart, idxEnd;
  UInt              degree;

public:
  ParametrizationTemporalPolynomial(Config &config);

  Bool setInterval(const Time &timeStart, const Time &timeEnd, Bool estimatePerArc);
  UInt parameterCount() const {return (degree+1)*(idxEnd-idxStart);}
  void factors(const Time &time, UInt startIndex, std::vector<UInt> &index, std::vector<Double> &factor) const;
  void parameterName(std::vector<ParameterName> &name) const;
};

/***********************************************/

#endif
