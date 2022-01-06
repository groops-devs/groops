/***********************************************/
/**
* @file instrumentResample.cpp
*
* @brief Resample data to given time series using polynomial prediction or least squares polynomial fit.
*
* @author Andreas Kvas
* @author Sebastian Strasser
* @author Torsten Mayer-Guerr
* @date 2019-10-07
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program resamples \file{instrument data}{instrument} to a given
\configClass{timeSeries}{timeSeriesType} using a resampling \config{method} out of the following:

\begin{itemize}
\item \config{polynomial}: Polynomial prediction using a moving polynomial of \config{polynomialDegree}.
  The optimal polynomial is chosen based on the centricity of the data points around the resampling
  point and the distance to all polynomial data points. All polynomial data points must be within
  \config{maxDataPointRange}. Resampling points within \config{maxExtrapolationDistance} of the
  polynomial will be extrapolated.

  \fig{!hb}{0.5}{instrumentResample_polynomial}{fig:instrumentResample_polynomial}{Example of polynomial prediction when resampling from 5 to 1 minute sampling}

\item \config{leastSquaresPolynomialFit}. A polynomial of \config{polynomialDegree} is estimated using
  all data points within \config{maxDataPointDistance} of the resampling point. This polynomial is then
  used to predict the resampling point. A resampling point will be extrapolated if there are only data
  points before/after as long as the closest one is within \config{maxExtrapolationDistance}.

  \fig{!hb}{0.5}{instrumentResample_leastSquares}{fig:instrumentResample_leastSquares}{Example of least squares polynomial fit when resampling from 5 to 1 minute sampling}
\end{itemize}

The elements \config{maxDataPointRange}, \config{maxDataPointDistance}, and \config{maxExtrapolationDistance}
are given in the unit of seconds. If negative values are used, the unit is relative to the median input sampling.

This program can also be used to reduce the sampling of an instrument file,
but a better way to reduce the sampling of noisy data with regular sampling
is to use a low pass filter first with \program{InstrumentFilter} and then thin
out the data with \program{InstrumentReduceSampling}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/polynomial.h"
#include "files/fileInstrument.h"
#include "classes/timeSeries/timeSeries.h"

/***** CLASS ***********************************/

/** @brief Resample data to given time series using polynomial prediction or least squares polynomial fit.
* @ingroup programsGroup */
class InstrumentResample
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(InstrumentResample, SINGLEPROCESS, "Resample data to given time series using polynomial prediction or least squares polynomial fit.", Instrument)

/***********************************************/

void InstrumentResample::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName      outName, inName;
    UInt          polynomialDegree;
    Double        maxDataPointRange, maxExtrapolationDistance;
    TimeSeriesPtr timeSeries;
    Bool          isLeastSquares = FALSE;
    std::string   choice;

    readConfig(config, "outputfileInstrument", outName, Config::MUSTSET,  "", "");
    readConfig(config, "inputfileInstrument",  inName,  Config::MUSTSET,  "", "");
    if(readConfigChoice(config, "method", choice, Config::MUSTSET,  "polynomial", "resampling method"))
    {
      if(readConfigChoiceElement(config, "polynomial", choice, "polynomial prediction"))
      {
        isLeastSquares = FALSE;
        readConfig(config, "polynomialDegree",         polynomialDegree,         Config::MUSTSET, "3",        "degree of the moving polynomial");
        readConfig(config, "maxDataPointRange",        maxDataPointRange,        Config::MUSTSET, "-(3+1.1)", "[seconds] all degree+1 data points must be within this range for a valid polynomial");
        readConfig(config, "maxExtrapolationDistance", maxExtrapolationDistance, Config::DEFAULT, "-1.1",     "[seconds] resampling points within this range of the polynomial will be extrapolated");
      }
      if(readConfigChoiceElement(config, "leastSquaresPolynomialFit", choice, "least squares polynomial fit"))
      {
        isLeastSquares = TRUE;
        readConfig(config, "polynomialDegree",         polynomialDegree,         Config::MUSTSET, "",  "degree of the estimated polynomial");
        readConfig(config, "maxDataPointDistance",     maxDataPointRange,        Config::MUSTSET, "",  "[seconds] all data points within this distance around the resampling point will be used");
        readConfig(config, "maxExtrapolationDistance", maxExtrapolationDistance, Config::DEFAULT, "0", "[seconds] resampling points within this range of the polynomial will be extrapolated");
      }
      endChoice(config);
    }
    readConfig(config, "timeSeries", timeSeries, Config::MUSTSET, "", "resampled points in time");
    if(isCreateSchema(config)) return;

    logStatus<<"read instrument data <"<<inName<<">"<<Log::endl;
    const Arc arc = InstrumentFile::read(inName);
    const std::vector<Time> timesNew = timeSeries->times();

    logStatus<<"resample data"<<Log::endl;
    Polynomial polynomial;
    polynomial.init(arc.times(), polynomialDegree, FALSE/*throwException*/, isLeastSquares, maxDataPointRange, maxExtrapolationDistance);
    const Matrix A = polynomial.interpolate(timesNew, arc.matrix().column(1, Epoch::dataCount(arc.getType(), TRUE)));

    Arc arcNew;
    Epoch *epoch = Epoch::create(arc.getType());
    for(UInt i=0; i<timesNew.size(); i++)
      if(!std::isnan(A(i,0)))
      {
        epoch->time  = timesNew.at(i);
        epoch->setData(A.row(i).trans());
        arcNew.push_back(*epoch);
    }
    delete epoch;

    logStatus<<"write instrument data to file <"<<outName<<">"<<Log::endl;
    InstrumentFile::write(outName, arcNew);
    Arc::printStatistics(arcNew);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
