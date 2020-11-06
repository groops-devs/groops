/***********************************************/
/**
* @file timeSeriesOrbitRevolutions.h
*
* @brief Read orbit file and create a time stamp for each ascending equator crossing.
* @see TimeSeries
*
* @author Norbert Zehentner
* @date 2017-05-16
*
*/
/***********************************************/

#ifndef __GROOPS_TIMESERIESORBITREVOLUTIONS__
#define __GROOPS_TIMESERIESORBITREVOLUTIONS__

// Latex documentation
#ifdef DOCSTRING_TimeSeries
static const char *docstringTimeSeriesOrbitRevolutions = R"(
\subsection{Revolution}
Reads an \file{orbit file}{instrument} and create a time stamp for each ascending equator crossing.
The time series can be restricted to the interval
starting from \config{timeStart} and before \config{timeEnd}.
)";
#endif

/***********************************************/

#include "classes/timeSeries/timeSeries.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Read orbit file and create a time stamp for each ascending equator crossing.
* @ingroup timeSeriesGroup
* @see TimeSeries */
class TimeSeriesOrbitRevolutions : public TimeSeriesBase
{
  FileName fileName;
  Time     timeStart, timeEnd;

public:
  TimeSeriesOrbitRevolutions(Config &config);

  std::vector<Time> times() const;
};

/***********************************************/

inline TimeSeriesOrbitRevolutions::TimeSeriesOrbitRevolutions(Config &config)
{
  try
  {
    timeEnd = date2time(9999, 1, 1);

    readConfig(config, "inputfileOrbit", fileName,  Config::MUSTSET,  "", "");
    readConfig(config, "timeStart",      timeStart, Config::OPTIONAL, "", "exclude peochs before this epoch");
    readConfig(config, "timeEnd",        timeEnd,   Config::OPTIONAL, "", "only epochs before this time are used");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline std::vector<Time> TimeSeriesOrbitRevolutions::times() const
{
  try
  {
    const OrbitArc orbit = InstrumentFile::read(fileName);

    std::vector<Time> times;
    for(UInt i=0; i<orbit.size()-1; i++)
      if((orbit.at(i).position.phi() < 0) && (orbit.at(i+1).position.phi() > 0))
      {
        const Double deltaPhi = Double(orbit.at(i+1).position.phi())-Double(orbit.at(i).position.phi());
        const Time   time     = orbit.at(i).time + (Double(orbit.at(i).position.phi())/deltaPhi) * (orbit.at(i+1).time-orbit.at(i).time);
        if(time.isInInterval(timeStart, timeEnd))
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
