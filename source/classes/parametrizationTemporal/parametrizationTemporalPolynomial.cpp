/***********************************************/
/**
* @file parametrizationTemporalPolynomial.cpp
*
* @brief Legendre Polynomials in time intervals.
* @see ParametrizationTemporal
*
* @author Torsten Mayer-Guerr
* @date 2015-05-31
*
*/
/***********************************************/

#include "base/import.h"
#include "base/legendrePolynomial.h"
#include "config/config.h"
#include "classes/timeSeries/timeSeries.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "parametrizationTemporalPolynomial.h"

/***********************************************/

ParametrizationTemporalPolynomial::ParametrizationTemporalPolynomial(Config &config)
{
  try
  {
    TimeSeriesPtr timeSeries;
    readConfig(config, "polynomialDegree", degree,     Config::MUSTSET,  "0", "polynomial degree");
    readConfig(config, "interval",         timeSeries, Config::DEFAULT,  "",  "intervals of polynomials");
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

Bool ParametrizationTemporalPolynomial::setInterval(const Time &timeStart, const Time &timeEnd, Bool estimatePerArc)
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

void ParametrizationTemporalPolynomial::parameterName(std::vector<ParameterName> &name) const
{
  try
  {
    for(UInt i=idxStart; i<idxEnd; i++)
    {
      name.push_back(ParameterName("", "", "", times.at(i), times.at(i+1))); // constant part
      for(UInt n=1; n<=degree; n++)
        name.push_back(ParameterName("", "", "legendrePolynomial.n"+n%"%i"s, times.at(i), times.at(i+1)));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ParametrizationTemporalPolynomial::factors(const Time &time, UInt startIndex, std::vector<UInt> &index, std::vector<Double> &factor) const
{
  try
  {
    if((time < times.at(idxStart)) || (time >= times.at(idxEnd)))
      return;

    // find index (interval)
    UInt idx = idxStart;
    while(time>=times.at(idx+1))
      idx++;

    const Double t  = 2*(time-times.at(idx)).mjd()/(times.at(idx+1)-times.at(idx)).mjd()-1; // t is [-1,+1]
    const Vector Pn = LegendrePolynomial::compute(t, degree);
    for(UInt n=0; n<=degree; n++)
    {
      index.push_back( (degree+1)*(idx-idxStart)+n+startIndex );
      factor.push_back( Pn(n) );
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
