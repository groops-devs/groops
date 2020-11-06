/***********************************************/
/**
* @file timeSeriesInstrumentArcIntervals.h
*
* @brief Reconstruct arc intervals from instrument file.
* @see TimeSeries
*
* @author Matthias Ellmer
* @date 2017-08-02
*
*/
/***********************************************/

#ifndef __GROOPS_TIMESERIESINSTRUMENTARCINTERVALS__
#define __GROOPS_TIMESERIESINSTRUMENTARCINTERVALS__

// Latex documentation
#ifdef DOCSTRING_TimeSeries
static const char *docstringTimeSeriesInstrumentArcIntervals = R"(
\subsection{InstrumentArcIntervals}
Reconstruct a time series from an \file{instrument file}{instrument}.
The time series is the first epoch of each arc plus one time step beyond the last
epoch of the last arc (using median sampling).
)";
#endif

/***********************************************/

#include "classes/timeSeries/timeSeries.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Read time series from a instrument file.
 * @ingroup timeSeriesGroup
 * @see TimeSeries */
class TimeSeriesInstrumentArcIntervals : public TimeSeriesBase
{
  FileName fileName;

public:
  TimeSeriesInstrumentArcIntervals(Config &config);

  std::vector<Time> times() const;
};

/***********************************************/

inline TimeSeriesInstrumentArcIntervals::TimeSeriesInstrumentArcIntervals(Config &config)
{
  try
  {
    readConfig(config, "inputfileInstrument", fileName,  Config::MUSTSET,  "", "Must be regular. Time series is first epoch of each arc plus one time step extrapolated from last epoch of last arc.");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline std::vector<Time> TimeSeriesInstrumentArcIntervals::times() const
{
  try
  {
    std::vector<Time> times;

    InstrumentFile file(fileName);
    Time finalTime;
    for(UInt arcNo=0; arcNo<file.arcCount(); arcNo++)
    {
      std::vector<Time> arcTimes = file.readArc(arcNo).times();
      if(arcTimes.size())
      {
        times.push_back(arcTimes.front());
        finalTime = arcTimes.back() + medianSampling(arcTimes);
      }
      else
        times.push_back(finalTime);
    }
    times.push_back(finalTime);

    return times;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif /* __GROOPS_TIMESERIESINSTRUMENTARCINTERVALS__ */
