/***********************************************/
/**
* @file timeSeriesEveryYear.h
*
* @brief Time series oriented on years.
* @see TimeSeries
*
* @author Torsten Mayer-Guerr
* @date 2017-11-06
*
*/
/***********************************************/

#ifndef __GROOPS_TIMESERIESYEVERYYEAR__
#define __GROOPS_TIMESERIESYEVERYYEAR__

// Latex documentation
#ifdef DOCSTRING_TimeSeries
static const char *docstringTimeSeriesEveryYear = R"(
\subsection{EveryYear}
Generates a time series with yearly sampling. The first point in time will be \config{timeStart} and the following
points are generated for each year at the same day and time in year.
The last generated point in time will be less or equal \config{timeEnd}.
)";
#endif

/***********************************************/

#include "classes/timeSeries/timeSeries.h"

/***** CLASS ***********************************/

/** @brief Time series oriented on years.
* @ingroup timeSeriesGroup
* @see TimeSeries */
class TimeSeriesEveryYear : public TimeSeriesBase
{
  Time timeStart, timeEnd;

public:
  TimeSeriesEveryYear(Config &config);

  std::vector<Time> times() const;
};

/***********************************************/

inline TimeSeriesEveryYear::TimeSeriesEveryYear(Config &config)
{
  try
  {
    readConfig(config, "timeStart", timeStart, Config::MUSTSET, "", "first point in time");
    readConfig(config, "timeEnd",   timeEnd,   Config::MUSTSET, "", "last point in time will be less or equal timeEnd");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline std::vector<Time> TimeSeriesEveryYear::times() const
{
  try
  {
    std::vector<Time> times = {timeStart};

    for(;;)
    {
      UInt   year, month, day, hour, min;
      Double sec;
      times.back().date(year, month, day, hour, min, sec);
      const Time time = date2time(year+1, month, day, hour, min, sec);
      if(time > timeEnd)
        break;
      times.push_back(time);
    }
    return times;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif /* __GROOPS__ */
