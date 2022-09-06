/***********************************************/
/**
* @file interpolatorTimeSeries.cpp
*
* @brief Interpolation of time series.
*
* @author Torsten Mayer-Guerr
* @date 2022-08-01
*
*/
/***********************************************/

#define DOCSTRING_InterpolatorTimeSeries

#include "base/import.h"
#include "config/configRegister.h"
#include "inputOutput/logging.h"
#include "classes/interpolatorTimeSeries/interpolatorTimeSeriesPolynomial.h"
#include "classes/interpolatorTimeSeries/interpolatorTimeSeriesLeastSquaresPolynomialFit.h"
#include "classes/interpolatorTimeSeries/interpolatorTimeSeriesFillGapsLeastSquaresPolynomialFit.h"
#include "classes/interpolatorTimeSeries/interpolatorTimeSeries.h"

/***********************************************/

GROOPS_REGISTER_CLASS(InterpolatorTimeSeries, "interpolatorTimeSeriesType",
                      InterpolatorTimeSeriesPolynomial,
                      InterpolatorTimeSeriesLeastSquaresPolynomialFit,
                      InterpolatorTimeSeriesFillGapsLeastSquaresPolynomialFit)

GROOPS_READCONFIG_CLASS(InterpolatorTimeSeries, "interpolatorTimeSeriesType")

/***********************************************/

InterpolatorTimeSeriesPtr InterpolatorTimeSeries::create(Config &config, const std::string &name)
{
  try
  {
    InterpolatorTimeSeriesPtr interpolatorTimeSeries;
    std::string choice;

    readConfigChoice(config, name, choice, Config::MUSTSET, "", "resampling method");
    if (readConfigChoiceElement(config, "polynomial",                        choice, "polynomial prediction"))
      interpolatorTimeSeries = InterpolatorTimeSeriesPtr(new InterpolatorTimeSeriesPolynomial(config));
    if (readConfigChoiceElement(config, "leastSquaresPolynomialFit",         choice, "least squares polynomial fit"))
      interpolatorTimeSeries = InterpolatorTimeSeriesPtr(new InterpolatorTimeSeriesLeastSquaresPolynomialFit(config));
    if (readConfigChoiceElement(config, "fillGapsLeastSquaresPolynomialFit", choice, "least squares polynomial fit"))
      interpolatorTimeSeries = InterpolatorTimeSeriesPtr(new InterpolatorTimeSeriesFillGapsLeastSquaresPolynomialFit(config));
    endChoice(config);

    return interpolatorTimeSeries;
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
