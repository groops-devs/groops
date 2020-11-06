/***********************************************/
/**
* @file timeSeriesUniformSampling.h
*
* @brief Time series with uniform sampling.
* @see TimeSeries
*
* @author Torsten Mayer-Guerr
* @date 2007-03-02
*
*/
/***********************************************/

#ifndef __GROOPS_TIMESERIESUNIFORMSAMPLING__
#define __GROOPS_TIMESERIESUNIFORMSAMPLING__

// Latex documentation
#ifdef DOCSTRING_TimeSeries
static const char *docstringTimeSeriesUniformSampling = R"(
\subsection{UniformSampling}\label{timeSeriesType:uniformSampling}
Generates a time series with uniform sampling. The first point in time will be \config{timeStart}.
The last generated point in time will be less or equal \config{timeEnd}.
The time step between generated points in time is given by \config{sampling}.
)";
#endif

/***********************************************/

#include "classes/timeSeries/timeSeries.h"

/***** CLASS ***********************************/

/** @brief Time series with uniform sampling.
* @ingroup timeSeriesGroup
* @see TimeSeries */
class TimeSeriesUniformSampling : public TimeSeriesBase
{
  Time timeStart, timeEnd, deltaTime;

public:
  TimeSeriesUniformSampling(Config &config);

  std::vector<Time> times() const;
};

/***********************************************/

inline TimeSeriesUniformSampling::TimeSeriesUniformSampling(Config &config)
{
  try
  {
    readConfig(config, "timeStart", timeStart, Config::MUSTSET, "", "first point in time");
    readConfig(config, "timeEnd",   timeEnd,   Config::MUSTSET, "", "last point in time will be less or equal timeEnd");
    readConfig(config, "sampling",  deltaTime, Config::MUSTSET, "", "time step between points in time");
    if(isCreateSchema(config)) return;

    if(deltaTime <= Time(0, 0.0))
      throw(Exception("Sampling must be strictly positive."));

  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline std::vector<Time> TimeSeriesUniformSampling::times() const
{
  try
  {
    std::vector<Time> times;
    for(UInt i=0; timeStart+i*deltaTime <= timeEnd; i++)
      times.push_back(timeStart+i*deltaTime);

    if(times.size()<1)
      throw(Exception("Points in time must be given in increasing order."));

    return times;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif /* __GROOPS__ */
