/***********************************************/
/**
* @file instrumentResample.cpp
*
* @brief Resample data to given time series using polynomial prediction or least squares polynomial fit.
*
* @author Andreas Kvas
* @author Sebastian Strasser
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

This program can also be used to reduce the sampling of an instrument file,
but a better way to reduce the sampling of noisy data with regular sampling
is to use a low pass filter first with \program{InstrumentFilter} and then thin
out the data with \program{InstrumentReduceSampling}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "classes/timeSeries/timeSeries.h"

/***** CLASS ***********************************/

/** @brief Resample data to given time series using polynomial prediction or least squares polynomial fit.
* @ingroup programsGroup */
class InstrumentResample
{
public:
  class ResamplerBase
  {
  public:
    virtual ~ResamplerBase() {}
    virtual Arc resample(Arc &arcOld, const std::vector<Time> &resamplingTimes) const = 0;
  };

  class ResamplerPolynomial : public ResamplerBase
  {
    private:
      Double maxDataPointRange, maxExtrapolationDistance;
      UInt polynomialDegree;

    public:
      ResamplerPolynomial(Config &config)
      {
        readConfig(config, "polynomialDegree",         polynomialDegree,         Config::MUSTSET,  "",  "degree of the moving polynomial");
        readConfig(config, "maxDataPointRange",        maxDataPointRange,        Config::MUSTSET,  "",  "[seconds] all degree+1 data points must be within this range for a valid polynomial");
        readConfig(config, "maxExtrapolationDistance", maxExtrapolationDistance, Config::DEFAULT,  "0", "[seconds] resampling points within this range of the polynomial will be extrapolated");
        if(isCreateSchema(config)) return;
      }
      Arc resample(Arc &arcOld, const std::vector<Time> &resamplingTimes) const override;
  };

  class ResamplerLeastSquares : public ResamplerBase
  {
    private:
      Double maxDataPointDistance, maxExtrapolationDistance;
      UInt polynomialDegree;

    public:
      ResamplerLeastSquares(Config &config)
      {
        readConfig(config, "polynomialDegree",         polynomialDegree,         Config::MUSTSET,  "",  "degree of the estimated polynomial");
        readConfig(config, "maxDataPointDistance",     maxDataPointDistance,     Config::MUSTSET,  "",  "[seconds] all data points within this distance around the resampling point will be used");
        readConfig(config, "maxExtrapolationDistance", maxExtrapolationDistance, Config::DEFAULT,  "0", "[seconds] resampling points within this range of the polynomial will be extrapolated");
        if(isCreateSchema(config)) return;
      }
      Arc resample(Arc &arcOld, const std::vector<Time> &resamplingTimes) const override;
  };

  class Resampler
  {
    private:
      std::shared_ptr<ResamplerBase> resampler;

    public:
      void initFromConfig(Config &config, const std::string &name, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
      {
        std::string choice;
        if(readConfigChoice(config, name, choice, mustSet, defaultValue, annotation))
        {
          if(readConfigChoiceElement(config, "polynomial", choice, "polynomial prediction"))
            resampler = std::shared_ptr<ResamplerBase>(new ResamplerPolynomial(config));
          if(readConfigChoiceElement(config, "leastSquaresPolynomialFit", choice, "least squares polynomial fit"))
            resampler = std::shared_ptr<ResamplerBase>(new ResamplerLeastSquares(config));
          endChoice(config);
        }
      }

      Arc resample(Arc &arcOld, const std::vector<Time> &resamplingTimes)
      {
        if(resampler == nullptr)
          return Arc();
        return resampler->resample(arcOld, resamplingTimes);
      }
  };

  static Vector predictByPolynomialFit(const std::vector<Time> &timesOld, const_MatrixSliceRef dataOld, const Time &resamplingTime, UInt polynomialDegree);
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(InstrumentResample, SINGLEPROCESS, "Resample data to given time series using polynomial prediction or least squares polynomial fit.", Instrument)

/***********************************************/

template<> Bool readConfig(Config &config, const std::string &name, InstrumentResample::Resampler &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  var.initFromConfig(config, name, mustSet, defaultValue, annotation);
  return TRUE;
}

/***********************************************/

void InstrumentResample::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName      outName, inName;
    TimeSeriesPtr timeSeries;
    Resampler     resampler;

    readConfig(config, "outputfileInstrument", outName,     Config::MUSTSET,  "", "");
    readConfig(config, "inputfileInstrument",  inName,      Config::MUSTSET,  "", "");
    readConfig(config, "method",               resampler,   Config::MUSTSET,  "polynomial", "resampling method");
    readConfig(config, "timeSeries",           timeSeries,  Config::MUSTSET,  "",           "resampled points in time");
    if(isCreateSchema(config)) return;

    // ======================================================

    logStatus<<"read instrument data <"<<inName<<">"<<Log::endl;
    Arc arc = InstrumentFile::read(inName);

    // create time series
    // ------------------
    logStatus<<"create time series"<<Log::endl;
    std::vector<Time> times = timeSeries->times();
    logInfo<<"  time start "<<times.at(0).dateTimeStr()<<Log::endl;
    logInfo<<"  time end   "<<times.back().dateTimeStr()<<Log::endl;

    // resample
    // -----------
    logStatus<<"resample data"<<Log::endl;
    Arc arcNew = resampler.resample(arc, times);

    // save file
    // ---------
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

Vector InstrumentResample::predictByPolynomialFit(const std::vector<Time> &timesOld, const_MatrixSliceRef dataOld, const Time &resamplingTime, UInt polynomialDegree)
{
  try
  {
    const UInt count = timesOld.size();

    // polynomial prediction or adjustment
    // -----------------------------------
    Matrix W(count, polynomialDegree+1);
    for(UInt k=0; k<count; k++)
    {
      const Double factor = (timesOld.at(k)-resamplingTime).seconds();
      W(k,0) = 1.0;
      for(UInt n=1; n<=polynomialDegree; n++)
        W(k,n) = factor * W(k,n-1);
    }

    // solve with QR - decomposition
    // -----------------------------
    Vector tau = QR_decomposition(W);
    Vector coeff(count);
    coeff(0) = 1.;
    triangularSolve(1., W.row(0, W.columns()).trans(), coeff.row(0, W.columns())); // R^(-T)*coeff
    QMult(W, tau, coeff); // coeff := Q*R^(-T)*coeff

    // apply prediction coefficients
    // -----------------------------
    return dataOld.trans() * coeff;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Arc InstrumentResample::ResamplerPolynomial::resample(Arc &arcOld, const std::vector<Time> &resamplingTimes) const
{
  try
  {
    const Matrix A = arcOld.matrix();
    const std::vector<Time> timesOld = arcOld.times();

    Arc arcNew;

    const Time dataPointRange = seconds2time(maxDataPointRange);
    const Time extrapolationDistance = seconds2time(maxExtrapolationDistance);

    auto searchStart = timesOld.begin();
    Single::forEach(resamplingTimes.size(), [&](UInt i)
    {
      const Time resamplingTime = resamplingTimes.at(i);

      // data point search space for potential polynomials
      searchStart    = std::lower_bound(searchStart, timesOld.end(), resamplingTime-dataPointRange-extrapolationDistance); // first epoch greater or equal than search interval
      auto searchEnd = std::upper_bound(searchStart, timesOld.end(), resamplingTime+dataPointRange+extrapolationDistance); // first epoch outside search interval

      // search for optimal polynomial by moving polynomial through search space from left to right
      auto   optimalStart      = searchEnd;
      UInt   optimalCentricity = MAX_UINT; // primary metric
      Double optimalDistance   = (maxDataPointRange+maxExtrapolationDistance)*(polynomialDegree+1); // secondary metric
      for(auto polyFirst = searchStart; polyFirst+polynomialDegree<searchEnd; polyFirst++)
      {
        const auto polyLast = polyFirst+polynomialDegree;

        if((*polyLast-*polyFirst) > dataPointRange ||         // polynomial data points not within maxDataPointRange
           resamplingTime < *polyFirst-extrapolationDistance ||  // outside maxExtrapolationDistance before polynomial
           resamplingTime > *polyLast+extrapolationDistance)  // outside maxExtrapolationDistance after polynomial
          continue;

        // distances from resampling point to all polynomial data points
        std::vector<Double> distances(polynomialDegree+1);
        for(UInt i=0; i<distances.size(); i++)
          distances.at(i) = (*(polyFirst+i)-resamplingTime).seconds();

        // primary metric for how centric the resampling point is within the polynomial data points (0 = centric, e.g. sum of {-1 -1 1 1})
        const UInt centricity = std::abs(std::accumulate(distances.begin(), distances.end(), 0, [](Int sum, Double distance){ return sum+(distance < 0. ? -1 : 1); }));
        if(centricity > optimalCentricity)
          break; // there won't be any better polynomial to the right of the current optimal anymore if centricity is increasing

        // secondary metric for how close the resampling point is to all polynomial data points
        const Double distance = std::fabs(std::accumulate(distances.begin(), distances.end(), 0.));

        if(centricity < optimalCentricity || (centricity == optimalCentricity && distance <= optimalDistance))
        {
          optimalStart      = polyFirst;
          optimalCentricity = centricity;
          optimalDistance   = distance;
        }
      }

      if(optimalStart == searchEnd)
        return; // no valid polynomial found

      // polynomial prediction
      Vector data = predictByPolynomialFit(std::vector<Time>(optimalStart, optimalStart+polynomialDegree+1),
                                           A.slice(std::distance(timesOld.begin(), optimalStart), 1, polynomialDegree+1, A.columns()-1),
                                           resamplingTime,
                                           polynomialDegree);

      searchStart = optimalStart; // limit search space on the left to the latest polynomial start data point

      Epoch *epoch = Epoch::create(arcOld.getType());
      epoch->time = resamplingTime;
      epoch->setData(data);
      arcNew.push_back(*epoch);
      delete epoch;
    });

    return arcNew;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}


/***********************************************/

Arc InstrumentResample::ResamplerLeastSquares::resample(Arc &arcOld, const std::vector<Time> &resamplingTimes) const
{
  try
  {
    const Matrix A = arcOld.matrix();
    const std::vector<Time> timesOld = arcOld.times();

    Arc arcNew;

    const Time dataPointDistance = seconds2time(maxDataPointDistance);
    const Time extrapolationDistance = seconds2time(maxExtrapolationDistance);

    auto searchStart = timesOld.begin();
    Single::forEach(resamplingTimes.size(), [&](UInt i)
    {
      const Time resamplingTime = resamplingTimes.at(i);

      searchStart    = std::lower_bound(searchStart, timesOld.end(), resamplingTime-dataPointDistance); // first epoch greater or equal than search interval
      auto searchEnd = std::upper_bound(searchStart, timesOld.end(), resamplingTime+dataPointDistance); // first epoch outside search interval

      std::vector<Time> deltaTime;
      std::transform(searchStart, searchEnd, std::back_inserter(deltaTime), [resamplingTime](const Time &t) { return t-resamplingTime; } );

      UInt count = deltaTime.size();
      if(count < (this->polynomialDegree+1)) // not enough points
        return;

      // check if we are extrapolating and if so if we are allowed to
      if( (deltaTime.back() < Time(0, 0) && (deltaTime.back()+extrapolationDistance) < Time(0, 0) ) || // all points are before resamplingTime and we are too far away or
          (Time(0, 0) < deltaTime.front() && (extrapolationDistance-deltaTime.front()) < Time(0, 0)) )         // all points are after resamplingTime and we are too far away
            return;

      // polynomial adjustment
      Vector data = predictByPolynomialFit(std::vector<Time>(searchStart, searchEnd),
                                           A.slice(std::distance(timesOld.begin(), searchStart), 1, count, A.columns()-1),
                                           resamplingTime,
                                           polynomialDegree);

      Epoch *epoch = Epoch::create(arcOld.getType());
      epoch->time = resamplingTime;
      epoch->setData(data);
      arcNew.push_back(*epoch);
      delete epoch;
    });

    return arcNew;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
