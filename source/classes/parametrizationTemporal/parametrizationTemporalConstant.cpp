/***********************************************/
/**
* @file parametrizationTemporalConstant.cpp
*
* @brief Constant in time.
* @see ParametrizationTemporal
*
* @author Torsten Mayer-Guerr
* @date 2011-12-28
*
*/
/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "classes/timeSeries/timeSeries.h"
#include "parametrizationTemporalConstant.h"

/***********************************************/

ParametrizationTemporalConstant::ParametrizationTemporalConstant(Config &config)
{
  try
  {
    TimeSeriesPtr timeSeries;
    readConfig(config, "interval",        timeSeries,      Config::DEFAULT, "", "");
    readConfig(config, "includeLastTime", includeLastTime, Config::DEFAULT, "0", "");
    if(isCreateSchema(config)) return;

    times      = timeSeries->times();
    isInterval = (times.size() != 0);
    if(!isInterval)
      times = {Time(), date2time(2500,1,1)};
    idxStart = 0;
    idxEnd   = times.size()-1;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool ParametrizationTemporalConstant::setInterval(const Time &timeStart, const Time &timeEnd, Bool estimatePerArc)
{
  try
  {
    const UInt idxStartOld = idxStart;
    const UInt idxEndOld   = idxEnd;

    if(estimatePerArc && !isInterval)
      times = {timeStart, timeEnd};

    idxStart = 0;
    while((idxStart+1<times.size()) && (timeStart>=times.at(idxStart+1)))
      idxStart++;
    idxEnd = idxStart;
    while((idxEnd<times.size()-1) && (timeEnd>times.at(idxEnd)))
      idxEnd++;

    return (idxStartOld != idxStart) || (idxEndOld != idxEnd);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ParametrizationTemporalConstant::parameterName(std::vector<ParameterName> &name) const
{
  try
  {
    for(UInt i=idxStart; i<idxEnd; i++)
      name.push_back(ParameterName("", "", "", times.at(i), times.at(i+1)));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ParametrizationTemporalConstant::factors(const Time &time, UInt startIndex, std::vector<UInt> &index, std::vector<Double> &factor) const
{
  try
  {
    if(!times.size())
    {
      index.push_back(startIndex);
      factor.push_back(1.);
      return;
    }

    if((time < times.at(idxStart)) || (time > times.at(idxEnd)) || (!includeLastTime && (time == times.at(idxEnd))))
      return;

    // find index (interval)
    UInt idx = idxStart;
    while((idx+1 < idxEnd) && (time >= times.at(idx+1)))
      idx++;

    index.push_back(idx-idxStart+startIndex);
    factor.push_back(1.);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
