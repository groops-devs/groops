/***********************************************/
/**
* @file timeSeriesUniformInterval.h
*
* @brief Time series from given interval count.
* @see TimeSeries
*
* @author Torsten Mayer-Guerr
* @date 2007-03-02
*
*/
/***********************************************/

#ifndef __GROOPS_TIMESERIESUNIFORMINTERVAL__
#define __GROOPS_TIMESERIESUNIFORMINTERVAL__

// Latex documentation
#ifdef DOCSTRING_TimeSeries
static const char *docstringTimeSeriesUniformInterval = R"(
\subsection{UniformInterval}
Generates a time series with uniform sampling between \config{timeStart} and \config{timeEnd}.
\config{intervallCount} gives the count of intervals. This class generates count+1 points in time
inclusive \config{timeStart} and \config{timeEnd}.
)";
#endif

/***********************************************/

#include "classes/timeSeries/timeSeries.h"

/***** CLASS ***********************************/

/** @brief Time series from given interval count.
* @ingroup timeSeriesGroup
* @see TimeSeries */
class TimeSeriesUniformInterval : public TimeSeriesBase
{
  Time timeStart, timeEnd;
  UInt count;

public:
  TimeSeriesUniformInterval(Config &config);

  std::vector<Time> times() const;
};

/***********************************************/

inline TimeSeriesUniformInterval::TimeSeriesUniformInterval(Config &config)
{
  try
  {
    readConfig(config, "timeStart",      timeStart, Config::MUSTSET, "", "1st point of the time series");
    readConfig(config, "timeEnd",        timeEnd,   Config::MUSTSET, "", "last point of the time series");
    readConfig(config, "intervalCount",  count,     Config::MUSTSET, "", "count of intervals, count+1 points in time will generated");
    if(isCreateSchema(config)) return;

    if(timeStart>=timeEnd)
      throw(Exception("Points in time must be given in increasing order."));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline std::vector<Time> TimeSeriesUniformInterval::times() const
{
  try
  {
    std::vector<Time> times(count+1);
    const Time deltaT = 1./count * (timeEnd-timeStart);
    for(UInt i=0; i<=count; i++)
      times.at(i) = timeStart + i * deltaT;
    return times;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif /* __GROOPS__ */
