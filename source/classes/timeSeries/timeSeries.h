/***********************************************/
/**
* @file timeSeries.h
*
* @brief Generates time series.
*
* @author Torsten Mayer-Guerr
* @date 2007-03-02
*
*/
/***********************************************/

#ifndef __GROOPS_TIMESERIES__
#define __GROOPS_TIMESERIES__

// Latex documentation
#ifdef DOCSTRING_TimeSeries
static const char *docstringTimeSeries = R"(
\section{TimeSeries}\label{timeSeriesType}
This class generates a series of points in time. The series is always sorted in ascending order.
Depending of the application the series is interpreted as list of points or as intervals between the points.

\fig{!hb}{0.4}{timeSeriesIntervals}{fig:timeSeriesIntervals}{List of points $t_i$ vs. intervals $T_i$.}
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"

/**
* @defgroup timeSeriesGroup TimeSeries
* @brief Generates time series.
* @ingroup classesGroup
* The interface is given by @ref TimeSeries. */
/// @{

/***** TYPES ***********************************/

class TimeSeries;
class TimeSeriesBase;
typedef std::shared_ptr<TimeSeries> TimeSeriesPtr;

/***** CLASS ***********************************/

/** @brief Generates time series.
* An instance of this class can be created with @ref readConfig. */
class TimeSeries
{
  std::vector<TimeSeriesBase*> base;

public:
  /// Constructor.
  TimeSeries(Config &config, const std::string &name);

  /// Destructor.
  ~TimeSeries();

  /** @brief Time series with increasing order. */
  std::vector<Time> times() const;

  /** @brief creates an derived instance of this class. */
  static TimeSeriesPtr create(Config &config, const std::string &name) {return TimeSeriesPtr(new TimeSeries(config, name));}
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class TimeSeries.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and a class without times is created.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] timeSeries Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates TimeSeries */
template<> Bool readConfig(Config &config, const std::string &name, TimeSeriesPtr &timeSeries, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***** CLASS ***********************************/

// Internal class
class TimeSeriesBase
{
public:
  virtual ~TimeSeriesBase() {}
  virtual std::vector<Time> times() const = 0;
};

/***********************************************/

#endif

