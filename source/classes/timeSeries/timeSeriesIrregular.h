/***********************************************/
/**
* @file timeSeriesIrregular.h
*
* @brief Time series given as explicit list.
* @see TimeSeries
*
* @author Torsten Mayer-Guerr
* @date 2007-03-02
*
*/
/***********************************************/

#ifndef __GROOPS_TIMESERIESIRREGULAR__
#define __GROOPS_TIMESERIESIRREGULAR__

// Latex documentation
#ifdef DOCSTRING_TimeSeries
static const char *docstringTimeSeriesIrregular = R"(
\subsection{Irregular}\label{timeSeriesType:irregular}
The points of the time series are given explicitly with \config{time}.
)";
#endif

/***********************************************/

#include "classes/timeSeries/timeSeries.h"

/***** CLASS ***********************************/

/** @brief Read time given as explicit list.
* @ingroup timeSeriesGroup
* @see TimeSeries */
class TimeSeriesIrregular : public TimeSeriesBase
{
  std::vector<Time> times_;

public:
  TimeSeriesIrregular(Config &config);

  std::vector<Time> times() const {return times_;}
};

/***********************************************/

inline TimeSeriesIrregular::TimeSeriesIrregular(Config &config)
{
  try
  {
    readConfig(config, "time", times_, Config::MUSTSET, "", "explicit list of points in time");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif /* __GROOPS__ */
