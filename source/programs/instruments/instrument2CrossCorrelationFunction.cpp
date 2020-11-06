/***********************************************/
/**
* @file instrument2CrossCorrelationFunction.cpp
*
* @brief Empirical computation of the correlation between two instrument files.
*
* @author Andreas Kvas
* @date 2018-01-01
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program computes the cross correlation between all corresponding data columns
in two \file{instrument files}{instrument}. The instrument files must be synchronized (\program{InstrumentSynchronize}).
The \configFile{outputfileCorrelation}{matrix} is a matrix with the first column containing the time lag followed by
cross-correlation function for each data column. The maximum lag is defined by the maximum arc length.

The correlation is based on the unbiased estimate of the cross-covariance between data columns $x$ and $y$,
\begin{equation}
  \sigma_{xy}(h) = \frac{1}{N}\sum_{k=1} x_{k+h} y_k,
\end{equation}
which is averaged over all arcs. From this estimate, the correlation for each lag is then computed via
\begin{equation}
  r_{xy}(h) = \frac{\sigma_{xy}(h)}{\sigma_x(0)\sigma_y(0)},
\end{equation}
which is the ratio between the biased estimates of the cross-covariance at lag $h$ and the auto-covariance of the individual data columns.

For instrument with data gaps, lag bins without any data are set to NAN.
)";

/***********************************************/

#include "programs/program.h"
#include "base/fourier.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Empirical computation of the correlation between two instrument files.
* @ingroup programsGroup */
class Instrument2CrossCorrelationFunction
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(Instrument2CrossCorrelationFunction, PARALLEL, "Empirical computation of the correlation between two instrument files.", Instrument, Statistics)
GROOPS_RENAMED_PROGRAM(InstrumentComputeCorrelation, Instrument2CrossCorrelationFunction, date2time(2020, 7, 7))

/***********************************************/

void Instrument2CrossCorrelationFunction::run(Config &config)
{
  try
  {
    FileName outputName;
    FileName inputName, inputNameReference;

    readConfig(config, "outputfileCorrelation",         outputName,         Config::MUSTSET, "", "column 1: time lag, column 2..n cross-correlation");
    readConfig(config, "inputfileInstrument",           inputName,          Config::MUSTSET,  "", "");
    readConfig(config, "inputfileInstrumentReference",  inputNameReference, Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    // check input
    // -----------
    InstrumentFile instrumentFile(inputName);
    InstrumentFile instrumentFileReference(inputNameReference);
    InstrumentFile::checkArcCount({instrumentFile, instrumentFileReference});

    const UInt arcCount  = instrumentFile.arcCount();
    const UInt dataCount = std::min(instrumentFile.dataCount(TRUE/*mustDefined*/), instrumentFileReference.dataCount(TRUE/*mustDefined*/));

    // determine arc length and data fields
    // ------------------------------------
    UInt   maxLag;
    Double sampling = 1.0;
    if(Parallel::isMaster())
    {
      std::vector<Time> times;
      Time maxArcLen;
      for(UInt arcNo = 0; arcNo<arcCount; arcNo++)
      {
        Arc arc = instrumentFile.readArc(arcNo);
        if(arc.size() == 0)
          continue;
        std::vector<Time> arcTimes = arc.times();
        maxArcLen = std::max(arcTimes.back() - arcTimes.front(), maxArcLen);
        times.insert(times.end(), arcTimes.begin(), arcTimes.end());
      }
      sampling = medianSampling(times).seconds();
      maxLag   = static_cast<UInt>(std::round(maxArcLen.seconds()/sampling)+1);
      logInfo<<"  maximum arc length: "<<maxLag<<" epochs"<<Log::endl;
      logInfo<<"  median sampling:    "<<sampling<<" seconds"<<Log::endl;
    }
    Parallel::broadCast(maxLag);
    Parallel::broadCast(sampling);

    // estimate the covariance for each arc, then reduce
    // -------------------------------------------------
    logStatus<<"Estimate cross-correlation for each arc"<<Log::endl;
    Vector count(maxLag*2-1);
    Vector autoCovarianceX(dataCount);
    Vector autoCovarianceY(dataCount);
    Matrix crossCovariance(maxLag*2-1, dataCount+1);
    Parallel::forEach(arcCount, [&](UInt arcNo)
    {
      const Arc arc    = instrumentFile.readArc(arcNo);
      const Arc arcRef = instrumentFileReference.readArc(arcNo);
      if(arc.size() == 0)
        return;

      const Matrix X = arc.matrix();
      const Matrix Y = arcRef.matrix();

      for(UInt i=0; i<dataCount; i++)
        autoCovarianceX(i) += quadsum(X.column(i+1));
      for(UInt i=0; i<dataCount; i++)
        autoCovarianceY(i) += quadsum(Y.column(i+1));

      std::vector<Time> times = arc.times();
      if(isRegular(times)) // fast version possible?
      {
        count(maxLag-1) += X.rows();
        for(UInt col=1; col<crossCovariance.columns(); col++)
          crossCovariance(maxLag-1, col) += inner(X.column(col), Y.column(col));
        for(UInt h=1; h<X.rows(); h++)
        {
          count(maxLag-1-h) += X.rows()-h;
          count(maxLag-1+h) += X.rows()-h;
          for(UInt col=1; col<crossCovariance.columns(); col++)
          {
            crossCovariance(maxLag-1-h, col) += inner(X.slice(0, col, X.rows()-h, 1), Y.slice(h, col, X.rows()-h, 1));
            crossCovariance(maxLag-1+h, col) += inner(X.slice(h, col, X.rows()-h, 1), Y.slice(0, col, X.rows()-h, 1));
          }
        }
      }
      else // general case
      {
        for(UInt i=0; i<X.rows(); i++)
          for(UInt k=0; k<Y.rows(); k++)
          {
            const UInt idx = static_cast<UInt>(std::round((times.at(i)-times.at(k)).seconds()/sampling)+maxLag-1);
            count(idx)++;
            for(UInt col=1; col<crossCovariance.columns(); col++)
              crossCovariance(idx, col) += X(i, col) * Y(k, col);
          }
      }
    });
    Parallel::reduceSum(count);
    Parallel::reduceSum(autoCovarianceX);
    Parallel::reduceSum(autoCovarianceY);
    Parallel::reduceSum(crossCovariance);

    if(Parallel::isMaster())
    {
      autoCovarianceX *= 1./count(maxLag-1);
      autoCovarianceY *= 1./count(maxLag-1);
      for(UInt i=0; i<count.rows(); i++)
        crossCovariance.row(i) *= 1./count(i);

      for(UInt i=0; i<dataCount; i++)
        crossCovariance.column(i+1) *= 1.0/std::sqrt(autoCovarianceX(i)*autoCovarianceY(i));

      for(UInt h=0; h<maxLag; h++)
      {
        crossCovariance(maxLag-h-1, 0) = -sampling*h;
        crossCovariance(h+maxLag-1, 0) =  sampling*h;
      }

      logStatus<<"write cross correlation to <"<<outputName<<">"<<Log::endl;
      writeFileMatrix(outputName, crossCovariance);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
