/***********************************************/
/**
* @file parametrizationTemporalFourier.h
*
* @brief Fourier expansion in time domain.
* @see ParametrizationTemporal
*
* @author Torsten Mayer-Guerr
* @date 2011-12-28
*
*/
/***********************************************/

#ifndef __GROOPS_PARAMETRIZATIONTEMPORALFOURIER__
#define __GROOPS_PARAMETRIZATIONTEMPORALFOURIER__

// Latex documentation
#ifdef DOCSTRING_ParametrizationTemporal
static const char *docstringParametrizationTemporalFourier = R"(
\subsection{Fourier}
A time variable function is given by a fourier expansion
\begin{equation}
f(x,t) = \sum_{m=1}^M f_m^c(\M x)\cos(2\pi m \tau) + f_m^s(\M x)\sin(2\pi m \tau)
\end{equation}
with the normalized time
\begin{equation}
\tau = \frac{t-t_A}{t_B-t_A},
\end{equation}
and $t_A$ is \config{timeStart}, $t_B$ is \config{timeEnd} in each \config{interval}
and $M$ is the \config{fourierDegree}.

The total parameter count is $2MN$, where $N$ is the number of intervals.
The parameters are sorted in following order: $f_1^c, f_1^s, f_2^c, \ldots$.
)";
#endif

/***********************************************/

#include "config/config.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"

/***** CLASS ***********************************/

/** @brief Fourier expansion.
* @ingroup parametrizationTemporalGroup
* @see ParametrizationTemporal */
class ParametrizationTemporalFourier : public ParametrizationTemporalBase
{
  std::vector<Time> times;
  Bool              isInterval;
  UInt              idxStart, idxEnd;
  UInt              order;

public:
  ParametrizationTemporalFourier(Config &config);

  Bool setInterval(const Time &timeStart, const Time &timeEnd, Bool estimatePerArc);
  UInt parameterCount() const  {return 2*order*(idxEnd-idxStart);}
  void parameterName(std::vector<ParameterName> &name) const;
  void factors(const Time &time, UInt startIndex, std::vector<UInt> &index, std::vector<Double> &factor) const;
};

/***********************************************/

#endif
