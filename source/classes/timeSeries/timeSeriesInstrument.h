/***********************************************/
/**
* @file timeSeriesInstrument.h
*
* @brief Read time series from a instrument file.
* @see TimeSeries
*
* @author Torsten Mayer-Guerr
* @date 2011-06-24
*
*/
/***********************************************/

#ifndef __GROOPS_TIMESERIESINSTRUMENT__
#define __GROOPS_TIMESERIESINSTRUMENT__

// Latex documentation
#ifdef DOCSTRING_TimeSeries
static const char *docstringTimeSeriesInstrument = R"(
\subsection{Instrument}\label{timeSeriesType:instrument}
Read a time series (epochs) from an \file{instrument file}{instrument}.
The time series can be restricted to the interval
starting from \config{timeStart} and before \config{timeEnd}.
)";
#endif

/***********************************************/

#include "classes/timeSeries/timeSeries.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Read time series from a instrument file.
* @ingroup timeSeriesGroup
* @see TimeSeries */
class TimeSeriesInstrument : public TimeSeriesBase
{
  FileName fileName;
  Time     timeStart, timeEnd;

public:
  TimeSeriesInstrument(Config &config);

  std::vector<Time> times() const;
};

/***********************************************/

inline TimeSeriesInstrument::TimeSeriesInstrument(Config &config)
{
  try
  {
    timeEnd = date2time(9999, 1, 1);

    readConfig(config, "inputfileInstrument", fileName,  Config::MUSTSET,  "", "");
    readConfig(config, "timeStart",           timeStart, Config::OPTIONAL, "", "exclude peochs before this epoch");
    readConfig(config, "timeEnd",             timeEnd,   Config::OPTIONAL, "", "only epochs before this time are used");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline std::vector<Time> TimeSeriesInstrument::times() const
{
  try
  {
    const Arc arc = InstrumentFile::read(fileName);

    std::vector<Time> times;
    for(UInt i=0; i<arc.size(); i++)
    {
      if(arc.at(i).time >= timeEnd)
        break;
      if(arc.at(i).time >= timeStart)
        times.push_back(arc.at(i).time);
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
