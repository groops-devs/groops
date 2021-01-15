/***********************************************/
/**
* @file instrument2SpectralCoherence.cpp
*
* @brief Empirical computation of the spectral coherence between two instrument files.
*
* @author Andreas Kvas
* @date 2018-01-01
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program computes the spectral coherence between two \file{instrument files}{instrument}.

The (magnitude-squared) coherence is defined as
\begin{equation}
  C_{xy}(f) = \frac{|P_{xy}(f)|^2}{P_{xx}(f)P_{yy}(f)}
\end{equation}
and is a measure in the range [0, 1] for the similarity of the signals $x$ and $y$ in frequency domain.
$P_{xy}$ is the cross-spectral density between $x$ and $y$ and $P_{xx}$, $P_{yy}$ are auto-spectral densities.
Auto- and cross-spectral densities are computed using Lomb's method (see \program{Instrument2PowerSpectralDensity} for details).

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

The \configFile{outputfileCoherence}{matrix} contains a matrix with the frequency vector as first column,
the coherence for each instrument channel is saved in the following columns.
)";

/***********************************************/

#include "programs/program.h"
#include "base/fourier.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Empirical computation of the spectral coherence between two instrument files.
* @ingroup programsGroup */
class Instrument2SpectralCoherence
{
  Matrix designMatrix(const Vector &t, Double f, Bool isNyquist = FALSE);
  std::vector<std::vector<std::complex<Double>>> leastSquaresFourier(const Vector &freqs, const_MatrixSliceRef arcMatrix, Bool countEven);

public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Instrument2SpectralCoherence, PARALLEL, "Empirical computation of the spectral coherence between two instrument files.", Instrument, Statistics)
GROOPS_RENAMED_PROGRAM(InstrumentComputeSpectralCoherence, Instrument2SpectralCoherence, date2time(2020, 7, 7))

/***********************************************/

void Instrument2SpectralCoherence::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName outputName;
    FileName inputName, inputNameReference;

    readConfig(config, "outputfileCoherence",           outputName,            Config::MUSTSET, "", "column 1: frequency, column 2-n coherence");
    readConfig(config, "inputfileInstrument",           inputName,             Config::MUSTSET,  "", "");
    readConfig(config, "inputfileInstrumentReference",  inputNameReference,    Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    // check input
    // -----------
    InstrumentFile instrumentFile(inputName);
    InstrumentFile instrumentFileReference(inputNameReference);
    InstrumentFile::checkArcCount({instrumentFile, instrumentFileReference});

    // determine arc length and data fields
    // ------------------------------------
    Vector freqs;
    UInt arcEpochCount, dataCount, arcCount;
    Double sampling = 1.0;
    if(Parallel::isMaster(comm))
    {
      arcCount = instrumentFile.arcCount();
      std::vector<Time> times;
      Time maxArcLen = seconds2time(0.0);

      dataCount = 0;
      for(UInt arcNo = 0; arcNo<arcCount; arcNo++)
      {
        Arc arc = instrumentFile.readArc(arcNo);
        if(arc.size() == 0)
          continue;
        auto arcTimes = arc.times();

        dataCount = std::max(dataCount, arc.at(0).data().rows());
        maxArcLen = std::max(arcTimes.back() - arcTimes.front(), maxArcLen);
        times.insert(times.end(), arcTimes.begin(), arcTimes.end());
      }
      sampling = medianSampling(times).seconds();

      arcEpochCount = static_cast<UInt>(std::round(maxArcLen.seconds()/sampling)+1);

      freqs = Fourier::frequencies(arcEpochCount, sampling);
      logInfo<<"  maximum arc length: "<<arcEpochCount<<" epochs"<<Log::endl;
      logInfo<<"  median sampling:    "<<sampling<<" seconds"<<Log::endl;
    }
    Parallel::broadCast(freqs,         0, comm);
    Parallel::broadCast(arcEpochCount, 0, comm);
    Parallel::broadCast(dataCount,     0, comm);
    Parallel::broadCast(arcCount,      0, comm);
    Bool countEven = (arcEpochCount%2) == 0; // flag that determines how to handle the nyquist frequency, see fourier.h for details

    // estimate the covariance matrix for each arc, then reduce
    // --------------------------------------------------------
    logStatus<<"Estimate spectral coherence for each arc"<<Log::endl;
    Matrix Gxx(freqs.rows(), dataCount);
    Matrix Gyy(freqs.rows(), dataCount);

    Matrix GxyReal(freqs.rows(), dataCount);
    Matrix GxyImag(freqs.rows(), dataCount);

    Parallel::forEach(arcCount, [&](UInt arcNo)
    {
      Arc arc = instrumentFile.readArc(arcNo);
      Arc arcRef = instrumentFileReference.readArc(arcNo);

      Matrix X = arc.matrix();
      Matrix Y = arcRef.matrix();

      std::vector<std::vector<std::complex<Double>>> F = leastSquaresFourier(freqs, X, countEven);
      std::vector<std::vector<std::complex<Double>>> G = leastSquaresFourier(freqs, Y, countEven);

      // accumulate estimates
      for(UInt k = 0; k<std::min(F.size(), G.size()); k++)
      {
        for(UInt n = 0; n<freqs.rows(); n++)
        {
          auto c = F.at(k).at(n)*std::conj(G.at(k).at(n)); // cross PSD
          GxyReal(n, k) += c.real();
          GxyImag(n, k) += c.imag();
          Gxx(n, k) += std::abs(F.at(k).at(n)*std::conj(F.at(k).at(n))); // auto PSDs
          Gyy(n, k) += std::abs(G.at(k).at(n)*std::conj(G.at(k).at(n)));
        }
      }
    }, comm);

    Parallel::reduceSum(GxyReal, 0, comm);
    Parallel::reduceSum(GxyImag, 0, comm);
    Parallel::reduceSum(Gxx, 0, comm);
    Parallel::reduceSum(Gyy, 0, comm);

    if(Parallel::isMaster(comm))
    {
      Matrix coherence(freqs.rows(), dataCount+1); // first row is frequency
      copy(freqs, coherence.column(0));

      for(UInt n = 0; n<freqs.rows(); n++)
        for(UInt k=0; k<dataCount; k++)
          coherence(n, k+1) = (GxyReal(n,k)*GxyReal(n,k) + GxyImag(n,k)*GxyImag(n,k))/(Gxx(n, k)*Gyy(n,k)); // C =|Gxy|^2/(Gxx*Gyy)

      logStatus<<"write coherence to <"<<outputName<<">"<<Log::endl;
      writeFileMatrix(outputName, coherence);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix Instrument2SpectralCoherence::designMatrix(const Vector &t, Double f, Bool isNyquist)
{
  Matrix A;
  if(isNyquist)
  {
    A = Matrix(t.rows(), 1);
    for(UInt i = 0; i<A.rows(); i++)
      A(i, 0) = (i%2) == 0 ? 1.0 : -1.0;
  }
  else
  {
    A = Matrix(t.rows(), 2);
    for(UInt i = 0; i<A.rows(); i++)
    {
      A(i, 0) = std::sin(2*PI*f*t(i));
      A(i, 1) = std::cos(2*PI*f*t(i));
    }
  }

  return A;
}

/***********************************************/

std::vector<std::vector<std::complex<Double>>> Instrument2SpectralCoherence::leastSquaresFourier(const Vector &freqs, const_MatrixSliceRef arcMatrix, Bool countEven)
{
  Vector t(arcMatrix.column(0));
  t-=t(0);
  t*=86400.0; // mjd -> seconds

  std::vector< std::vector< std::complex<Double> > > F(arcMatrix.columns()-1);

  for(UInt i = 0; i<F.size(); i++) // zero frequency: mean
    F.at(i).push_back(std::complex<Double>(mean(arcMatrix.column(i+1)), 0.0));

  UInt loopCount = countEven ? freqs.size()-2 : freqs.size()-1;
  for(UInt k = 0; k<loopCount; k++)
  {
    Matrix A = designMatrix(t, freqs[k+1]);
    Matrix l = Matrix(arcMatrix.column(1, arcMatrix.columns()-1));

    Matrix x_hat = leastSquares(A, l);
    for(UInt i = 0; i<arcMatrix.columns()-1; i++)
      F.at(i).push_back(std::complex<Double>(0.5*x_hat(1, i), -0.5*x_hat(0, i)));
  }
  if(countEven) // special case nyquist frequency
  {
    Matrix A = designMatrix(t, 0.5, TRUE);
    Matrix l = Matrix(arcMatrix.column(1, arcMatrix.columns()-1));

    Matrix x_hat = leastSquares(A, l);
    for(UInt i = 0; i<arcMatrix.columns()-1; i++)
      F.at(i).push_back(std::complex<Double>(x_hat(0, i), 0.0));
  }

  return F;
}

/***********************************************/
