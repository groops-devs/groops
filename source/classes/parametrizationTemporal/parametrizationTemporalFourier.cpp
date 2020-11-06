/***********************************************/
/**
* @file parametrizationTemporalFourier.cpp
*
* @brief Fourier expansion.
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
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "parametrizationTemporalFourier.h"

/***********************************************/

ParametrizationTemporalFourier::ParametrizationTemporalFourier(Config &config)
{
  try
  {
    TimeSeriesPtr timeSeries;

    readConfig(config, "fourierDegree", order,      Config::MUSTSET,  "", "");
    readConfig(config, "interval",      timeSeries, Config::DEFAULT,  "", "intervals of fourier series");
    if(isCreateSchema(config)) return;

    times      = timeSeries->times();
    isInterval = (times.size() != 0);
    if(!isInterval)
      times = {Time(), date2time(2500,1,1)};
    idxStart   = 0;
    idxEnd = times.size()-1;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ParametrizationTemporalFourier::setInterval(const Time &timeStart, const Time &timeEnd, Bool estimatePerArc)
{
  try
  {
    if(estimatePerArc && !isInterval)
      times = {timeStart, timeEnd};

    idxStart = 0;
    while((idxStart+1<times.size()) && (timeStart>=times.at(idxStart+1)))
      idxStart++;
    idxEnd = idxStart;
    while((idxEnd<times.size()-1) && (timeEnd>times.at(idxEnd)))
      idxEnd++;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ParametrizationTemporalFourier::parameterName(std::vector<ParameterName> &name) const
{
  try
  {
    for(UInt i=idxStart; i<idxEnd; i++)
      for(UInt m=1; m<=order; m++)
      {
        name.push_back(ParameterName("", "", "fourier.cos("+m%"%i"s+"*x)", times.at(i), times.at(i+1)));
        name.push_back(ParameterName("", "", "fourier.sin("+m%"%i"s+"*x)", times.at(i), times.at(i+1)));
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ParametrizationTemporalFourier::factors(const Time &time, UInt startIndex, std::vector<UInt> &index, std::vector<Double> &factor) const
{
  try
  {
    if((time<times.at(idxStart))||(time>=times.at(idxEnd)))
      return;

    // find index (interval)
    UInt idx = idxStart;
    while(time>=times.at(idx+1))
      idx++;

    const Double angle = 2*PI*(time-times.at(idx)).mjd()/(times.at(idx+1)-times.at(idx)).mjd();
    for(UInt m=1; m<=order; m++)
    {
      index.push_back(2*order*(idx-idxStart)+2*m-2+startIndex);
      index.push_back(2*order*(idx-idxStart)+2*m-1+startIndex);
      factor.push_back(cos(m*angle));
      factor.push_back(sin(m*angle));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
