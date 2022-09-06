/***********************************************/
/**
* @file interpolatorTimeSeriesPolynomial.h
*
* @brief Interpolation using a moving polynomial.
*
* @author Torsten Mayer-Guerr
* @date 2022-08-01
*
*/
/***********************************************/

#ifndef __GROOPS_INTERPOLATORTIMESERIESPOLYNOMIAL__
#define __GROOPS_INTERPOLATORTIMESERIESPOLYNOMIAL__

// Latex documentation
#ifdef DOCSTRING_InterpolatorTimeSeries
static const char *docstringInterpolatorTimeSeriesPolynomial = R"(
\subsection{Polynomial}
Polynomial prediction using a moving polynomial of \config{polynomialDegree}.
The optimal polynomial is chosen based on the centricity of the data points around the resampling
point and the distance to all polynomial data points. All polynomial data points must be within
\config{maxDataPointRange}. Resampling points within \config{maxExtrapolationDistance} of the
polynomial will be extrapolated. The elements \config{maxDataPointRange} and \config{maxExtrapolationDistance}
are given in the unit of seconds. If negative values are used, the unit is relative to the median input sampling.

\fig{!hb}{0.5}{instrumentResample_polynomial}{fig:instrumentResample_polynomial}{Example of polynomial prediction when resampling from 5 to 1 minute sampling}
)";
#endif

/***********************************************/

#include "base/polynomial.h"
#include "classes/interpolatorTimeSeries/interpolatorTimeSeries.h"

/***** CLASS ***********************************/

/** @brief Interpolation using a moving polynomial.
* @ingroup interpolatorTimeSeriesGroup
* @see InterpolatorTimeSeries */
class InterpolatorTimeSeriesPolynomial : public InterpolatorTimeSeries
{
  Polynomial polynomial;
  UInt       polynomialDegree;
  Double     maxDataPointRange, maxExtrapolationDistance;

public:
  InterpolatorTimeSeriesPolynomial(Config &config);

  void init(const std::vector<Time> &times, Bool throwException) override;
  Matrix interpolate(const std::vector<Time> &timesNew, const_MatrixSliceRef A, UInt rowsPerEpoch) const override;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline InterpolatorTimeSeriesPolynomial::InterpolatorTimeSeriesPolynomial(Config &config)
{
  try
  {
    readConfig(config, "polynomialDegree",         polynomialDegree,         Config::MUSTSET, "3",        "degree of the moving polynomial");
    readConfig(config, "maxDataPointRange",        maxDataPointRange,        Config::MUSTSET, "-(3+1.1)", "[seconds] all degree+1 data points must be within this range for a valid polynomial");
    readConfig(config, "maxExtrapolationDistance", maxExtrapolationDistance, Config::DEFAULT, "-1.1",     "[seconds] resampling points within this range of the polynomial will be extrapolated");
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void InterpolatorTimeSeriesPolynomial::init(const std::vector<Time> &times, Bool throwException)
{
  try
  {
    polynomial.init(times, polynomialDegree, throwException, FALSE/*isLeastSquares*/, maxDataPointRange, maxExtrapolationDistance);
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Matrix InterpolatorTimeSeriesPolynomial::interpolate(const std::vector<Time> &timesNew, const_MatrixSliceRef A, UInt rowsPerEpoch) const
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
