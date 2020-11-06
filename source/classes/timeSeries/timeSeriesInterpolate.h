/***********************************************/
/**
* @file timeSeriesInterpolate.h
*
* @brief Interpolates between times of a series.
* @see TimeSeries
*
* @author Torsten Mayer-Guerr
* @date 2017-11-05
*
*/
/***********************************************/

#ifndef __GROOPS_TIMESERIESINTERPOLATE__
#define __GROOPS_TIMESERIESINTERPOLATE__

// Latex documentation
#ifdef DOCSTRING_TimeSeries
static const char *docstringTimeSeriesInterpolate = R"(
\subsection{Interpolate}
Interpolates \config{nodeInterpolation} count points between
the given \configClass{timeSeries}{timeSeriesType} uniformly.
)";
#endif

/***********************************************/

#include "classes/timeSeries/timeSeries.h"

/***** CLASS ***********************************/

/** @brief Interpolates between times of a series
* @ingroup timeSeriesGroup
* @see TimeSeries */
class TimeSeriesInterpolate : public TimeSeriesBase
{
  TimeSeriesPtr timeSeries;
  UInt          interpCount;

public:
  TimeSeriesInterpolate(Config &config);

  std::vector<Time> times() const;
};

/***********************************************/

inline TimeSeriesInterpolate::TimeSeriesInterpolate(Config &config)
{
  try
  {
    readConfig(config, "timeSeries",        timeSeries,  Config::MUSTSET, "",  "time series to be created");
    readConfig(config, "nodeInterpolation", interpCount, Config::MUSTSET, "1", "interpolates count points in each time interval given by the time series");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline std::vector<Time> TimeSeriesInterpolate::times() const
{
  try
  {
    std::vector<Time> times = timeSeries->times();

    const UInt timeCount = times.size();
    for(UInt i=0; i<timeCount-1; i++)
      for(UInt k=1; k<=interpCount; k++)
        times.push_back(times.at(i) + (k/(interpCount+1.))*(times.at(i+1)-times.at(i)));
    std::sort(times.begin(), times.end());

    return times;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif /* __GROOPS__ */
