/***********************************************/
/**
* @file parametrizationTemporalOscillation.cpp
*
* @brief Oscillations.
* @see ParametrizationTemporal
*
* @author Torsten Mayer-Guerr
* @date 2015-06-07
*
*/
/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "parametrizationTemporalOscillation.h"

/***********************************************/

ParametrizationTemporalOscillation::ParametrizationTemporalOscillation(Config &config)
{
  try
  {
    readConfig(config, "period", timePeriod, Config::MUSTSET, "365.25",     "[day]");
    readConfig(config, "time0",  time0,      Config::MUSTSET, STRING_J2000, "reference time");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ParametrizationTemporalOscillation::parameterName(std::vector<ParameterName> &name) const
{
  try
  {
    for(UInt i=0; i<timePeriod.size(); i++)
    {
      const std::string str = "(2*pi/"+timePeriod.at(i).mjd()%"%g"s+"*(t-"+time0.dateTimeStr()+"))";
      name.push_back(ParameterName("", "", "oscillation.cos"+str));
      name.push_back(ParameterName("", "", "oscillation.sin"+str));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ParametrizationTemporalOscillation::factors(const Time &time, UInt startIndex, std::vector<UInt> &index, std::vector<Double> &factor) const
{
  try
  {
    for(UInt i=0; i<timePeriod.size(); i++)
    {
      const Double omega = 2*PI/timePeriod.at(i).mjd()*(time-time0).mjd();
      index.push_back(2*i+0+startIndex);
      index.push_back(2*i+1+startIndex);
      factor.push_back(cos(omega));
      factor.push_back(sin(omega));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
