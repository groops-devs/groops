/***********************************************/
/**
* @file instrument2PowerSpectralDensity.cpp
*
* @brief Compute PSD from instrument files.
*
* @author Andreas Kvas
* @date 2016-02-02
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program computes the power spectral density (PSD) for all data fields in an \file{instrument file}{instrument}.
The PSD is computed using Lomb's method. For each arc and each frequency $f$, a sinusoid is fit to the data
\begin{equation}
  l_i = a \cos(2\pi f t_i) + b \sin(2\pi f t_i) + e_i
\end{equation}

The PSD for this frequency is then computed by forming the square sum of adjusted observations:
\begin{equation}
  P(f) = \sum_i \hat{l}^2_i.
\end{equation}

The resulting PSD is the average over all arcs. For regularly sampled time series,
this method yields the same results as FFT based PSD estimates.

A regular frequency grid based on the longest arc and the median sampling is computed.
The maximum number of epochs per arc is determined by
\begin{equation}
  N = \frac{t_{\text{end}} - t_{\text{start}}}{\Delta t_{\text{median}} } + 1,
\end{equation}
the Nyquist frequency is given by
\begin{equation}
  f_{\text{nyq}} = \frac{1}{2\Delta t_{\text{median}}}.
\end{equation}

If it is suspected that \configFile{inputfileInstrument}{instrument} contains secular variations,
the input should be detrended using \program{InstrumentDetrend}.

See also \program{Instrument2CovarianceFunctionVCE},
\program{CovarianceFunction2PowerSpectralDensity}, \program{PowerSpectralDensity2CovarianceFunction}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/fourier.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Compute PSD from instrument files.
* @ingroup programsGroup */
class Instrument2PowerSpectralDensity
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(Instrument2PowerSpectralDensity, PARALLEL, "Compute PSD from instrument files.", Instrument, Covariance)
GROOPS_RENAMED_PROGRAM(InstrumentComputePSD, Instrument2PowerSpectralDensity, date2time(2018, 7, 18))
GROOPS_RENAMED_PROGRAM(InstrumentComputePsd, Instrument2PowerSpectralDensity, date2time(2020, 7, 7))

/***********************************************/

void Instrument2PowerSpectralDensity::run(Config &config)
{
  try
  {
    FileName fileNameInstrument, fileNamePSD;

    readConfig(config, "outputfilePSD",        fileNamePSD,         Config::MUSTSET,  "", "estimated PSD: column 0: frequency vector, column 1-(n-1): PSD estimate for each channel");
    readConfig(config, "inputfileInstrument",  fileNameInstrument,  Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    logStatus<<"read instrument data"<<Log::endl;
    InstrumentFile instrumentFile(fileNameInstrument);
    const UInt arcCount  = instrumentFile.arcCount();
    const UInt dataCount = instrumentFile.dataCount(TRUE/*mustDefined*/);

    Vector freqs;
    UInt   arcEpochCount;
    Double sampling = 1.0;
    if(Parallel::isMaster())
    {
      std::vector<Time> times;
      Time maxArcLen;
      for(UInt arcNo=0; arcNo<arcCount; arcNo++)
      {
        Arc arc = instrumentFile.readArc(arcNo);
        if(arc.size() == 0)
          continue;
        std::vector<Time> arcTimes = arc.times();
        maxArcLen = std::max(arcTimes.back() - arcTimes.front(), maxArcLen);
        times.insert(times.end(), arcTimes.begin(), arcTimes.end());
      }
      sampling      = medianSampling(times).seconds();
      arcEpochCount = static_cast<UInt>(std::round(maxArcLen.seconds()/sampling)+1);
      logInfo<<"  maximum arc length: "<<arcEpochCount<<" epochs"<<Log::endl;
      logInfo<<"  median sampling:    "<<sampling<<" seconds"<<Log::endl;

      freqs = Fourier::frequencies(arcEpochCount, sampling);
    }
    Parallel::broadCast(freqs);
    Parallel::broadCast(arcEpochCount);

    logStatus<<"compute PSD"<<Log::endl;
    Matrix PSD(freqs.rows(), dataCount+1);
    Parallel::forEach(arcCount, [&](UInt arcNo)
    {
      Arc arc = instrumentFile.readArc(arcNo);
      Matrix data = arc.matrix();

      // time vector
      Vector t(arc.size());
      for(UInt i=0; i<t.rows(); i++)
        t(i) = (arc.at(i).time-arc.at(0).time).seconds();

      // square sum of observations
      Vector lPl(data.columns()-1);
      for(UInt i=0; i<lPl.rows(); i++)
        lPl(i) = quadsum(data.column(i+1));

      // estimate the power of each frequency
      for(UInt k=0; k<freqs.size(); k++)
      {
        Matrix l = data.column(1, data.columns()-1);
        Matrix A(l.rows(), 2);
        const Double f = 2*PI*freqs.at(k);
        for(UInt i=0; i<A.rows(); i++)
        {
          A(i, 0) = std::cos(f*t(i));
          A(i, 1) = std::sin(f*t(i));
        }
        if((freqs.at(k) == 0) || (std::fabs(freqs.at(k)*sampling-0.5) < 1e-5)) // zero or nyquist freq?
          A = A.column(0);

        reduceLeastSquaresFit(A, l); // l = e_hat
        for(UInt i=0; i<l.columns(); i++)
          PSD(k, i+1) += lPl(i) - quadsum(l.column(i));
      }
    });
    Parallel::reduceSum(PSD);

    if(Parallel::isMaster())
    {
      PSD *= sampling/arcCount; // PSD unit: input^2/Hz
      copy(freqs, PSD.column(0));
      logStatus<<"write PSD to file <"<<fileNamePSD<<">"<<Log::endl;
      writeFileMatrix(fileNamePSD, PSD);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
