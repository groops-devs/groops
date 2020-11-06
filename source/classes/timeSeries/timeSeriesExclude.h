/***********************************************/
/**
* @file timeSeriesExclude.h
*
* @brief Exclude times from a given time series.
* @see TimeSeries
*
* @author Torsten Mayer-Guerr
* @date 2017-11-05
*
*/
/***********************************************/

#ifndef __GROOPS_TIMESERIESEXCLUDE__
#define __GROOPS_TIMESERIESEXCLUDE__

// Latex documentation
#ifdef DOCSTRING_TimeSeries
static const char *docstringTimeSeriesExclude = R"(
\subsection{Exclude}
In a first step a \configClass{timeSeries}{timeSeriesType} is generated.
In a second step all times are removed which are in range before or after \config{excludeMargin} seconds
of the times given by \configClass{timeSeriesExclude}{timeSeriesType}.
)";
#endif

/***********************************************/

#include "classes/timeSeries/timeSeries.h"

/***** CLASS ***********************************/

/** @brief Exclude times from a given time series.
* @ingroup timeSeriesGroup
* @see TimeSeries */
class TimeSeriesExclude : public TimeSeriesBase
{
  TimeSeriesPtr timeSeries, timeSeriesExcl;
  Double        margin;

public:
  TimeSeriesExclude(Config &config);

  std::vector<Time> times() const;
};

/***********************************************/

inline TimeSeriesExclude::TimeSeriesExclude(Config &config)
{
  try
  {
    readConfig(config, "timeSeries",        timeSeries,     Config::MUSTSET,  "",     "time series to be created");
    readConfig(config, "timeSeriesExclude", timeSeriesExcl, Config::MUSTSET,  "",     "exclude this time points from time series (within margin)");
    readConfig(config, "excludeMargin",     margin,         Config::DEFAULT,  "1e-5", "on both sides [seconds]");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline std::vector<Time> TimeSeriesExclude::times() const
{
  try
  {
    std::vector<Time> times;
    std::vector<Time> timesOld  = timeSeries->times();
    std::vector<Time> timesExcl = timeSeriesExcl->times();

    UInt idxTime = 0;
    for(UInt i=0; i<timesOld.size(); i++)
    {
      while((idxTime < timesExcl.size()) && ((timesExcl.at(idxTime)-timesOld.at(i)).seconds() < -margin))
        idxTime++;
      if((idxTime >= timesExcl.size()) || ((timesExcl.at(idxTime)-timesOld.at(i)).seconds() > +margin))
        times.push_back(timesOld.at(i));
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
