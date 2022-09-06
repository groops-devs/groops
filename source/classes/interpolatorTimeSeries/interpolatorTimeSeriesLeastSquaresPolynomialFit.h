/***********************************************/
/**
* @file interpolatorTimeSeriesLeastSquaresPolynomialFit.h
*
* @brief Interpolation using a least squares polynomial fit.
*
* @author Torsten Mayer-Guerr
* @date 2022-08-01
*
*/
/***********************************************/

#ifndef __GROOPS_INTERPOLATORTIMESERIESLEASTSQUARESPOLYNOMIALFIT__
#define __GROOPS_INTERPOLATORTIMESERIESLEASTSQUARESPOLYNOMIALFIT__

// Latex documentation
#ifdef DOCSTRING_InterpolatorTimeSeries
static const char *docstringInterpolatorTimeSeriesLeastSquaresPolynomialFit = R"(
\subsection{Least squares polynomial fit}
A polynomial of \config{polynomialDegree} is estimated using all data points within
\config{maxDataPointDistance} of the resampling point. This polynomial is then used
to predict the resampling point. A resampling point will be extrapolated if there are
only data points before/after as long as the closest one is within \config{maxExtrapolationDistance}.
The elements \config{maxDataPointDistance} and \config{maxExtrapolationDistance} are given
in the unit of seconds. If negative values are used, the unit is relative to the median input sampling.

\fig{!hb}{0.5}{instrumentResample_leastSquares}{fig:instrumentResample_leastSquares}{Example of least squares polynomial fit when resampling from 5 to 1 minute sampling}
)";
#endif

/***********************************************/

#include "base/polynomial.h"
#include "classes/interpolatorTimeSeries/interpolatorTimeSeries.h"

/***** CLASS ***********************************/

/** @brief Interpolation using a least squares polynomial fit.
* @ingroup interpolatorTimeSeriesGroup
* @see InterpolatorTimeSeries */
class InterpolatorTimeSeriesLeastSquaresPolynomialFit : public InterpolatorTimeSeries
{
  Polynomial polynomial;
  UInt       polynomialDegree;
  Double     maxDataPointDistance, maxExtrapolationDistance;

public:
  InterpolatorTimeSeriesLeastSquaresPolynomialFit(Config &config);

  void init(const std::vector<Time> &times, Bool throwException) override;
  Matrix interpolate(const std::vector<Time> &timesNew, const_MatrixSliceRef A, UInt rowsPerEpoch) const override;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline InterpolatorTimeSeriesLeastSquaresPolynomialFit::InterpolatorTimeSeriesLeastSquaresPolynomialFit(Config &config)
{
  try
  {
    readConfig(config, "polynomialDegree",         polynomialDegree,         Config::MUSTSET, "",  "degree of the estimated polynomial");
    readConfig(config, "maxDataPointDistance",     maxDataPointDistance,     Config::MUSTSET, "",  "[seconds] all data points within this distance around the resampling point will be used");
    readConfig(config, "maxExtrapolationDistance", maxExtrapolationDistance, Config::DEFAULT, "0", "[seconds] resampling points within this range of the polynomial will be extrapolated");
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void InterpolatorTimeSeriesLeastSquaresPolynomialFit::init(const std::vector<Time> &times, Bool throwException)
{
  try
  {
    polynomial.init(times, polynomialDegree, throwException, TRUE/*isLeastSquares*/, maxDataPointDistance, maxExtrapolationDistance);
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Matrix InterpolatorTimeSeriesLeastSquaresPolynomialFit::interpolate(const std::vector<Time> &timesNew, const_MatrixSliceRef A, UInt rowsPerEpoch) const
{
  try
  {
    return polynomial.interpolate(timesNew, A, rowsPerEpoch);
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
