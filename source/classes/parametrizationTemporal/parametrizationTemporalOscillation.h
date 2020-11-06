/***********************************************/
/**
* @file parametrizationTemporalOscillation.h
*
* @brief Oscillations.
* @see ParametrizationTemporal
*
* @author Torsten Mayer-Guerr
* @date 2015-06-07
*
*/
/***********************************************/

#ifndef __GROOPS_PARAMETRIZATIONTEMPORALOSCILLATION__
#define __GROOPS_PARAMETRIZATIONTEMPORALOSCILLATION__

// Latex documentation
#ifdef DOCSTRING_ParametrizationTemporal
static const char *docstringParametrizationTemporalOscillation = R"(
\subsection{Oscillation}
A time variable function is given by a oscillation
\begin{equation}
f(x,t) = f^c(\M x)\cos(\omega_i(t)) + f^s(\M x)\sin(\omega_i(t))
\end{equation}
with $\omega_i=\frac{2\pi}{T_i}(t-t_0)$,
$t_0$ is \config{timeStart} and $T$ is \config{timePeriod} in days.
)";
#endif

/***********************************************/

#include "config/config.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"

/***** CLASS ***********************************/

/** @brief Oscillations.
* @ingroup parametrizationTemporalGroup
* @see ParametrizationTemporal */
class ParametrizationTemporalOscillation : public ParametrizationTemporalBase
{
  Time               time0;
  std::vector<Time>  timePeriod;

public:
  ParametrizationTemporalOscillation(Config &config);

  void setInterval(const Time &/*timeStart*/, const Time &/*timeEnd*/, Bool /*estimatePerArc*/) {}
  UInt parameterCount() const  {return 2*timePeriod.size();}
  void parameterName(std::vector<ParameterName> &name) const;
  void factors(const Time &time, UInt startIndex, std::vector<UInt> &index, std::vector<Double> &factor) const;
};

/***********************************************/

#endif
