/***********************************************/
/**
* @file timeSeriesMonthly.h
*
* @brief Time series oriented on months.
* @see TimeSeries
*
* @author Torsten Mayer-Guerr
* @date 2007-03-02
*
*/
/***********************************************/

#ifndef __GROOPS_TIMESERIESYMONTHLY__
#define __GROOPS_TIMESERIESYMONTHLY__

// Latex documentation
#ifdef DOCSTRING_TimeSeries
static const char *docstringTimeSeriesMonthly = R"(
\subsection{Monthly}
If \config{monthMiddle} is set, time points are generated at mid of each month inclusively
the \config{monthStart} in \config{yearStart} and \config{monthEnd} in \config{yearEnd}.
Otherwise times are given at the first of each month and a time point after the last month.
)";
#endif

/***********************************************/

#include "classes/timeSeries/timeSeries.h"

/***** CLASS ***********************************/

/** @brief Time series oriented on months.
* @ingroup timeSeriesGroup
* @see TimeSeries */
class TimeSeriesMonthly : public TimeSeriesBase
{
  UInt monthStart, yearStart;
  UInt monthEnd,   yearEnd;
  Bool monthMiddle;

public:
  TimeSeriesMonthly(Config &config);

  std::vector<Time> times() const;
};

/***********************************************/

inline TimeSeriesMonthly::TimeSeriesMonthly(Config &config)
{
  try
  {
    readConfig(config, "monthStart",     monthStart,  Config::MUSTSET,  "", "");
    readConfig(config, "yearStart",      yearStart,   Config::MUSTSET,  "", "");
    readConfig(config, "monthEnd",       monthEnd,    Config::MUSTSET,  "", "");
    readConfig(config, "yearEnd",        yearEnd,     Config::MUSTSET,  "", "");
    readConfig(config, "useMonthMiddle", monthMiddle, Config::DEFAULT,  "0", "time points are mid of months, otherwise the 1st of each month + a time point behind the last month");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline std::vector<Time> TimeSeriesMonthly::times() const
{
  try
  {
    std::vector<Time> times;
    for(UInt year=yearStart; year<=yearEnd; year++)
    {
      const UInt monthStart2 = (year==yearStart) ? monthStart : 1;
      const UInt monthEnd2   = (year==yearEnd)   ? monthEnd   : 12;
      for(UInt month=monthStart2; month<=monthEnd2; month++)
        if(monthMiddle)
          times.push_back( 0.5 * (date2time(year, month, 1) + date2time(year, month+1, 1)) );
        else
          times.push_back( date2time(year, month, 1) );
    }
    if(!monthMiddle)
      times.push_back(date2time(yearEnd, monthEnd+1, 1));
    return times;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif /* __GROOPS__ */
