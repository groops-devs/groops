/***********************************************/
/**
* @file parametrizationTemporalSplines.cpp
*
* @brief Splines expansion.
* @see ParametrizationTemporal
*
* @author Torsten Mayer-Guerr
* @date 2011-12-28
*
*/
/***********************************************/

#include "base/import.h"
#include "base/basisSplines.h"
#include "config/config.h"
#include "inputOutput/logging.h"
#include "classes/timeSeries/timeSeries.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "parametrizationTemporalSplines.h"

/***********************************************/

ParametrizationTemporalSplines::ParametrizationTemporalSplines(Config &config)
{
  try
  {
    TimeSeriesPtr timeSeriesPtr, timesIntervalPtr;

    readConfig(config, "degree",          degree,           Config::MUSTSET, "",  "degree of splines");
    readConfig(config, "timeSeries",      timeSeriesPtr,    Config::MUSTSET, "",  "nodal points in time domain");
    readConfig(config, "intervals",       timesIntervalPtr, Config::DEFAULT, "",  "");
    readConfig(config, "includeLastTime", includeLastTime,  Config::DEFAULT, "0", "");
    if(isCreateSchema(config)) return;

    times = timeSeriesPtr->times();

    // init time intervals
    // -------------------
    std::vector<Time> timesInterval = timesIntervalPtr->times();
    isInterval = (timesInterval.size() != 0);
    if(!isInterval)
      timesInterval = {times.at(0), times.back()};
    computeIntervals(timesInterval);
    idxStart = 0;
    idxEnd   = idEpochStart.size();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ParametrizationTemporalSplines::computeIntervals(const std::vector<Time> &timesInterval)
{
  try
  {
    idEpochStart.clear();
    idEpochEnd.clear();
    for(UInt idInterval=0; idInterval<timesInterval.size()-1; idInterval++)
    {
      UInt start = 0;
      while((start+1<times.size()) && (times.at(start+1)<=timesInterval.at(idInterval)))
        start++;
      UInt count = 0;
      while((start+count<times.size()) && (times.at(start+count)<=timesInterval.at(idInterval+1)))
        count++;
      if(count==0)
        continue;

      for(UInt n=0; n<degree; n++)
        idEpochStart.push_back(start);
      for(UInt i=0; i<count-1; i++)
      {
        idEpochStart.push_back(start+i);
        idEpochEnd.push_back(start+i+1);
      }
      for(UInt n=0; n<degree; n++)
        idEpochEnd.push_back(start+count-1);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool ParametrizationTemporalSplines::setInterval(const Time &timeStart, const Time &timeEnd, Bool estimatePerArc)
{
  try
  {
    const UInt idxStartOld = idxStart;
    const UInt idxEndOld   = idxEnd;

    if((timeEnd<=times.at(0)) || (timeStart>=times.back()))
    {
      idxStart = idxEnd = 0; // time interval outside spline nodal points
      return (idxStartOld != idxStart) || (idxEndOld != idxEnd);
    }

    if(estimatePerArc && !isInterval)
      computeIntervals({timeStart, timeEnd});

    idxStart = 0;
    while((idxStart<idEpochEnd.size()) && (times.at(idEpochEnd.at(idxStart))<=timeStart))
      idxStart++;

    idxEnd = idxStart;
    while((idxEnd<idEpochStart.size()) && (times.at(idEpochStart.at(idxEnd))<timeEnd))
      idxEnd++;

    return (idxStartOld != idxStart) || (idxEndOld != idxEnd);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ParametrizationTemporalSplines::parameterName(std::vector<ParameterName> &name) const
{
  try
  {
    for(UInt i=idxStart; i<idxEnd; i++)
      name.push_back(ParameterName("", "", "spline.n"+degree%"%i"s, times.at(idEpochStart.at(i)), times.at(idEpochEnd.at(i))));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ParametrizationTemporalSplines::factors(const Time &time, UInt startIndex, std::vector<UInt> &index, std::vector<Double> &factor) const
{
  try
  {
    if(idxEnd == idxStart)
      return;
    if((time < times.at(idEpochStart.at(idxStart))) || (time > times.at(idEpochEnd.at(idxEnd-1))) || (!includeLastTime && (time == times.at(idEpochEnd.at(idxEnd-1)))))
      return;

    // find first used parameter
    UInt idx = idxStart;
    while((idx+1+degree < idxEnd) && (time >= times.at(idEpochEnd.at(idx))))
      idx++;

    const UInt   idEpoch = idEpochEnd.at(idx)-1;
    const Double t = (time-times.at(idEpoch)).mjd()/(times.at(idEpoch+1)-times.at(idEpoch)).mjd();
    const Vector f = BasisSplines::compute(t, degree);
    for(UInt i=0; i<=degree; i++)
      if(std::fabs(f(i)) > 1e-8)
      {
        index.push_back( i+idx-idxStart+startIndex );
        factor.push_back( f(i) );
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
