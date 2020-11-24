/***********************************************/
/**
* @file parametrizationTemporalSplines.h
*
* @brief Splines in time domain.
* @see ParametrizationTemporal
*
* @author Torsten Mayer-Guerr
* @date 2011-12-28
*
*/
/***********************************************/

#ifndef __GROOPS_PARAMETRIZATIONTEMPORALSPLINES__
#define __GROOPS_PARAMETRIZATIONTEMPORALSPLINES__

// Latex documentation
#ifdef DOCSTRING_ParametrizationTemporal
static const char *docstringParametrizationTemporalSplines = R"(
\subsection{Splines}\label{parametrizationTemporalType:splines}
A time variable function is given by
\begin{equation}
f(x,t) =  \sum_i f_i(x)\Psi_i(t),
\end{equation}
with the (spatial) coefficients $f_i(x)$ as parameters and the temporal basis functions~$\Psi_i(t)$.
Basis splines are defined as polynomials of degree~$n$ in intervals between nodal points in time $t_i$,
for details see~\reference{basis splines}{fundamentals.basisSplines}.

The parameters are ordered timewise. First all parameters of $f_{i=1}(x)$ then
$f_{i=2}(x)$ and so on. The total parameter count in each \config{interval} is $N=N_t+d-1$,
where $N_t$ is the count of time points from \configClass{timeSeries}{timeSeriesType} in each interval and $d$
is the \config{degree}.
)";
#endif

/***********************************************/

#include "classes/parametrizationTemporal/parametrizationTemporal.h"

/***** CLASS ***********************************/

/** @brief Splines in time domain.
* @ingroup parametrizationTemporalGroup
* @see ParametrizationTemporal */
class ParametrizationTemporalSplines : public ParametrizationTemporalBase
{
  std::vector<Time> times;
  std::vector<UInt> idEpochStart, idEpochEnd;
  Bool              isInterval;
  UInt              idxStart, idxEnd;
  UInt              degree;

  void computeIntervals(const std::vector<Time> &timesInterval);

public:
  ParametrizationTemporalSplines(Config &config);

  Bool setInterval(const Time &timeStart, const Time &timeEnd, Bool estimatePerArc);
  UInt parameterCount() const {return idxEnd-idxStart;}
  void factors(const Time &time, UInt startIndex, std::vector<UInt> &index, std::vector<Double> &factor) const;
  void parameterName(std::vector<ParameterName> &name) const;
};

/***********************************************/

#endif
