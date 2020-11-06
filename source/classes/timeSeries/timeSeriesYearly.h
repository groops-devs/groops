/***********************************************/
/**
* @file timeSeriesYearly.h
*
* @brief Time series oriented on years.
* @see TimeSeries
*
* @author Torsten Mayer-Guerr
* @date 2016-06-08
*
*/
/***********************************************/

#ifndef __GROOPS_TIMESERIESYEARLY__
#define __GROOPS_TIMESERIESYEARLY__

// Latex documentation
#ifdef DOCSTRING_TimeSeries
static const char *docstringTimeSeriesYearly = R"(
\subsection{Yearly}
If \config{yearMiddle} is set, time points are generated at mid of each year inclusively \config{yearStart}
and \config{yearEnd}. Otherwise times are given at the first of each year and a time point after the last year.
)";
#endif

/***********************************************/

#include "classes/timeSeries/timeSeries.h"

/***** CLASS ***********************************/

/** @brief Time series oriented on years.
* @ingroup timeSeriesGroup
* @see TimeSeries */
class TimeSeriesYearly : public TimeSeriesBase
{
  UInt yearStart, yearEnd;
  Bool yearMiddle;

public:
  TimeSeriesYearly(Config &config);

  std::vector<Time> times() const;
};

/***********************************************/

inline TimeSeriesYearly::TimeSeriesYearly(Config &config)
{
  try
  {
    readConfig(config, "yearStart",         yearStart,   Config::MUSTSET,  "", "");
    readConfig(config, "yearEnd",           yearEnd,     Config::MUSTSET,  "", "");
    readConfig(config, "useYearMiddle",     yearMiddle,  Config::DEFAULT,  "0", "time points are mid of years, otherwise the 1st of each year + a time point behind the last year");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline std::vector<Time> TimeSeriesYearly::times() const
{
  try
  {
    std::vector<Time> times;
    for(UInt year=yearStart; year<=yearEnd; year++)
    {
      if(yearMiddle)
        times.push_back( 0.5 * (date2time(year, 1, 1) + date2time(year+1, 1, 1)) );
      else
        times.push_back( date2time(year, 1, 1) );
    }
    if(!yearMiddle)
      times.push_back(date2time(yearEnd+1, 1, 1));
    return times;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif /* __GROOPS__ */
