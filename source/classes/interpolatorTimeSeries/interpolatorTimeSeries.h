/***********************************************/
/**
* @file interpolatorTimeSeries.h
*
* @brief Interpolation of time series data.
*
* @author Torsten Mayer-Guerr
* @date 2022-08-01
*
*/
/***********************************************/

#ifndef __GROOPS_INTERPOLATORTIMESERIES__
#define __GROOPS_INTERPOLATORTIMESERIES__

// Latex documentation
#ifdef DOCSTRING_InterpolatorTimeSeries
static const char *docstringInterpolatorTimeSeries = R"(
\section{InterpolatorTimeSeries}\label{interpolatorTimeSeriesType}
This class resamples data of a times series to new poins in time.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"

/**
* @defgroup interpolatorTimeSeriesGroup InterpolatorTimeSeries
* @brief Interpolation of time series data.
* @ingroup classesGroup
* The interface is given by @ref InterpolatorTimeSeries.
* An Instance can be created by @ref readConfig. */
/// @{

/***** TYPES ***********************************/

class InterpolatorTimeSeries;
typedef std::shared_ptr<InterpolatorTimeSeries> InterpolatorTimeSeriesPtr;

/***** CLASS ***********************************/

/** @brief Interpolation of time series data.
* An Instance of this class can be created by @ref readConfig. */
class InterpolatorTimeSeries
{
public:
  /// Destructor.
  virtual ~InterpolatorTimeSeries() {}

  /** @brief Initialize the interpolatorTimeSeries.
  * @param times epochs of the input data.
  * @param throwException otherwise the non-interpolated values filled with NaN. */
  virtual void init(const std::vector<Time> &times, Bool throwException) = 0;

  /** @brief Interpolate a matrix to new epochs.
  * @param timesNew output epochs of the returned matrix.
  * @param A input data of time series
  * @param rowsPerEpoch e.g. for @a A with positions (x,y,z) per epoch in separated rows.
  * @return Interpolated matrix with timesNew.size()*rowsPerEpoch rows. */
  virtual Matrix interpolate(const std::vector<Time> &timesNew, const_MatrixSliceRef A, UInt rowsPerEpoch=1) const = 0;

  /** @brief creates a derived instance of this class. */
  static InterpolatorTimeSeriesPtr create(Config &config, const std::string &name);
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class InterpolatorTimeSeries.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a interpolatorTimeSeries is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] interpolatorTimeSeries Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates InterpolatorTimeSeries */
template<> Bool readConfig(Config &config, const std::string &name, InterpolatorTimeSeriesPtr &interpolatorTimeSeries, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

#endif
