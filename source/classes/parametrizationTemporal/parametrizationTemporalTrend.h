/***********************************************/
/**
* @file parametrizationTemporalTrend.h
*
* @brief Linear trend in time.
* @see ParametrizationTemporal
*
* @author Torsten Mayer-Guerr
* @date 2011-12-28
*
*/
/***********************************************/

#ifndef __GROOPS_PARAMETRIZATIONTEMPORALTREND__
#define __GROOPS_PARAMETRIZATIONTEMPORALTREND__

// Latex documentation
#ifdef DOCSTRING_ParametrizationTemporal
static const char *docstringParametrizationTemporalTrend = R"(
\subsection{Trend}\label{parametrizationTemporalType:trend}
A time variable function is given by a linear trend
\begin{equation}
f(x,t) = \frac{1}{T}(t-t_0) \cdot f_t(x),
\end{equation}
with $t_0$ is \config{timeStart} and $T$ is \config{timeStep} in days.
A constant term is not included and must added separately.
)";
#endif

/***********************************************/

#include "config/config.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"

/***** CLASS ***********************************/

/** @brief Linear trend in time.
* @ingroup parametrizationTemporalGroup
* @see ParametrizationTemporal */
class ParametrizationTemporalTrend : public ParametrizationTemporalBase
{
  Time time0, timeStep;

public:
  ParametrizationTemporalTrend(Config &config);

  Bool setInterval(const Time &/*timeStart*/, const Time &/*timeEnd*/, Bool /*estimatePerArc*/) {return FALSE;}
  UInt parameterCount() const {return 1;}
  void parameterName(std::vector<ParameterName> &name) const;
  void factors(const Time &time, UInt startIndex, std::vector<UInt> &index, std::vector<Double> &factor) const;
};

/***********************************************/

inline ParametrizationTemporalTrend::ParametrizationTemporalTrend(Config &config)
{
  try
  {
    readConfig(config, "timeStart", time0,    Config::MUSTSET,  STRING_J2000, "reference time");
    readConfig(config, "timeStep",  timeStep, Config::MUSTSET,  "365.25",     "");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void ParametrizationTemporalTrend::parameterName(std::vector<ParameterName> &name) const
{
  try
  {
    name.push_back(ParameterName("", "", "trend."+timeStep.mjd()%"%g"s+"*(t-"+time0.dateTimeStr()+")"));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void ParametrizationTemporalTrend::factors(const Time &time, UInt startIndex, std::vector<UInt> &index, std::vector<Double> &factor) const
{
  index.push_back(startIndex);
  factor.push_back((time-time0).mjd()/timeStep.mjd());
}

/***********************************************/

#endif
