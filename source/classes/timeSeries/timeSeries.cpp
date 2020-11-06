/***********************************************/
/**
* @file timeSeries.cpp
*
* @brief Generates time series.
*
* @author Torsten Mayer-Guerr
* @date 2007-03-02
*
*/
/***********************************************/

#define DOCSTRING_TimeSeries

#include "base/import.h"
#include "config/configRegister.h"
#include "classes/timeSeries/timeSeries.h"
#include "classes/timeSeries/timeSeriesUniformSampling.h"
#include "classes/timeSeries/timeSeriesUniformInterval.h"
#include "classes/timeSeries/timeSeriesIrregular.h"
#include "classes/timeSeries/timeSeriesMonthly.h"
#include "classes/timeSeries/timeSeriesYearly.h"
#include "classes/timeSeries/timeSeriesEveryMonth.h"
#include "classes/timeSeries/timeSeriesEveryYear.h"
#include "classes/timeSeries/timeSeriesInstrument.h"
#include "classes/timeSeries/timeSeriesInstrumentArcIntervals.h"
#include "classes/timeSeries/timeSeriesOrbitRevolutions.h"
#include "classes/timeSeries/timeSeriesExclude.h"
#include "classes/timeSeries/timeSeriesConditional.h"
#include "classes/timeSeries/timeSeriesInterpolate.h"

/***********************************************/

GROOPS_REGISTER_CLASS(TimeSeries, "timeSeriesType",
                      TimeSeriesUniformSampling,
                      TimeSeriesUniformInterval,
                      TimeSeriesIrregular,
                      TimeSeriesMonthly,
                      TimeSeriesYearly,
                      TimeSeriesEveryMonth,
                      TimeSeriesEveryYear,
                      TimeSeriesInstrument,
                      TimeSeriesInstrumentArcIntervals,
                      TimeSeriesOrbitRevolutions,
                      TimeSeriesExclude,
                      TimeSeriesConditional,
                      TimeSeriesInterpolate)

GROOPS_READCONFIG_UNBOUNDED_CLASS(TimeSeries, "timeSeriesType")

/***********************************************/

TimeSeries::TimeSeries(Config &config, const std::string &name)
{
  try
  {
    std::string type;
    while(readConfigChoice(config, name, type, Config::OPTIONAL, "", "Generates time series"))
    {
      std::vector<Time> times;
      if(readConfigChoiceElement(config, "uniformSampling",        type, "uniform series with given sampling"))
        base.push_back(new TimeSeriesUniformSampling(config));
      if(readConfigChoiceElement(config, "uniformInterval",        type, "uniform series with given interval count"))
        base.push_back(new TimeSeriesUniformInterval(config));
      if(readConfigChoiceElement(config, "irregular",              type, "every single point in time is given"))
        base.push_back(new TimeSeriesIrregular(config));
      if(readConfigChoiceElement(config, "monthly",                type, "time series based on months"))
        base.push_back(new TimeSeriesMonthly(config));
      if(readConfigChoiceElement(config, "yearly",                 type, "time series based on years"))
        base.push_back(new TimeSeriesYearly(config));
      if(readConfigChoiceElement(config, "everyMonth",             type, "time series based on months"))
        base.push_back(new TimeSeriesEveryMonth(config));
      if(readConfigChoiceElement(config, "everyYear",              type, "time series based on years"))
        base.push_back(new TimeSeriesEveryYear(config));
      if(readConfigChoiceElement(config, "instrument",             type, "read time series from an instrument file"))
        base.push_back(new TimeSeriesInstrument(config));
      if(readConfigChoiceElement(config, "instrumentArcIntervals", type, "try to reproduce arc interval time series from an instrument file"))
        base.push_back(new TimeSeriesInstrumentArcIntervals(config));
      if(readConfigChoiceElement(config, "revolutions",            type, "create time series from orbit file (revolutions)"))
        base.push_back(new TimeSeriesOrbitRevolutions(config));
      if(readConfigChoiceElement(config, "exclude",                type, "exclude times from a given time series"))
        base.push_back(new TimeSeriesExclude(config));
      if(readConfigChoiceElement(config, "conditional",            type, "time series based depending on conditions"))
        base.push_back(new TimeSeriesConditional(config));
      if(readConfigChoiceElement(config, "interpolate",            type, "interpolate between created times"))
        base.push_back(new TimeSeriesInterpolate(config));
      endChoice(config);
      if(isCreateSchema(config))
        return;
    };
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

TimeSeries::~TimeSeries()
{
  for(UInt i=0; i<base.size(); i++)
    delete base.at(i);
}

/***********************************************/

std::vector<Time> TimeSeries::times() const
{
  try
  {
    std::vector<Time> times;
    for(UInt i=0; i<base.size(); i++)
    {
      auto times2 = base.at(i)->times();
      times.insert(times.end(), times2.begin(), times2.end());
    }
    std::sort(times.begin(), times.end());
    return times;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
