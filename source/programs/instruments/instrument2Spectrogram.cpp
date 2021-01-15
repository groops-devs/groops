/***********************************************/
/**
* @file instrument2Spectrogram.cpp
*
* @brief compute spectrogram using Short Time Fourier Transform applied to an instrument file.
*
* @author Saniya Behzadpour
* @date 2018-07-05
*/

/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program applies the Short Time Fourier Transform (STFT) to selected data columns
of \configFile{inputfileInstrument}{instrument} and computes the spectrogram.
The STFT is computed at centered \configClass{timeSeries}{timeSeriesType} with
an (possible overlapping) rectangular window with \config{windowLength} seconds.
Data gaps are zero padded within the window.

The \configFile{outputfileSpectrogram}{matrix} is a matrix with each row the time (MJD),
the frequency $[Hz]$, and the amplitudes $[unit/\sqrt{Hz}]$ for the selected data columns.
It can be plotted with \program{PlotGraph}.

\fig{!hb}{0.8}{instrument2Spectrogram}{fig:instrument2Spectrogram}{GRACE range-rate residuals of one month (window of 6 hours).}
)";

/***********************************************/

#include "programs/program.h"
#include "base/fourier.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"
#include "classes/timeSeries/timeSeries.h"

/***** CLASS ***********************************/

/** @brief compute spectrogram using Short Time Fourier Transform applied to an instrument file.
* @ingroup programsGroup */
class Instrument2Spectrogram
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Instrument2Spectrogram, PARALLEL, "spectrogram using Short Time Fourier Transform", Instrument, Statistics)
GROOPS_RENAMED_PROGRAM(InstrumentComputeSpectrogram, Instrument2Spectrogram, date2time(2020, 7, 7))

/***********************************************/

void Instrument2Spectrogram::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName      fileNameOut, fileNameIn;
    TimeSeriesPtr timeSeries;
    Double        windowLength;
    UInt          startData, countData = MAX_UINT;

    readConfig(config, "outputfileSpectrogram", fileNameOut,  Config::MUSTSET, "",  "mjd, freq, ampl0, ampl1, ...");
    readConfig(config, "inputfileInstrument",   fileNameIn,   Config::MUSTSET, "",  "");
    readConfig(config, "timeSeries",            timeSeries,   Config::MUSTSET, "",  "center of SFFT window");
    readConfig(config, "windowLength",          windowLength, Config::MUSTSET, "",  "[seconds]");
    readConfig(config, "startDataFields",       startData,    Config::DEFAULT, "0", "start");
    readConfig(config, "countDataFields",       countData,    Config::OPTIONAL, "", "number of data fields (default: all)");
    if(isCreateSchema(config)) return;

    logStatus<<"read instrument data <"<<fileNameIn<<">"<<Log::endl;
    const Arc    arc         = InstrumentFile::read(fileNameIn);
    const Matrix A           = arc.matrix();
    const Double sampling    = medianSampling(arc.times()).seconds();
    const UInt   windowCount = static_cast<UInt>(std::round(windowLength/sampling));
    const Vector freqs       = Fourier::frequencies(windowCount, sampling);
    const std::vector<Time> times    = arc.times();
    const std::vector<Time> timesBin = timeSeries->times();
    logInfo<<" sampling = "<<sampling<<" seconds"<<Log::endl;
    logInfo<<" max. epoch in window = "<<windowCount<<Log::endl;

    countData = std::min(countData, A.columns()-1-startData);

    logStatus<<"computing spectrogram"<<Log::endl;
    Matrix spectrogram(timesBin.size()*freqs.rows(), 2+countData, NAN);
    auto iter = times.begin();
    Parallel::forEach(timesBin.size(), [&](UInt i)
    {
      const UInt idx = i*freqs.rows();
      copy(Vector(freqs.rows(), timesBin.at(i).mjd()), spectrogram.slice(idx, 0, freqs.rows(), 1)); // time
      copy(freqs, spectrogram.slice(idx, 1, freqs.rows(), 1)); // frequencies

      // find data in interval [timesBin-windowLength/2, timesBin+windowLength/2]
      iter = std::lower_bound(iter, times.end(), timesBin.at(i)-seconds2time(windowLength/2));
      const UInt row  = std::distance(times.begin(), iter);
      const UInt rows = std::distance(iter, std::upper_bound(iter, times.end(), timesBin.at(i)+seconds2time(windowLength/2)));
      if(rows == 0)
        return;

      // data with zero padding
      Matrix signal(windowCount, countData);
      copy(A.slice(row, 1+startData, std::min(rows, windowCount), countData), signal.row(0, std::min(rows, windowCount)));

      // compute fft
      for(UInt k=0; k<countData; k++)
      {
        std::vector<std::complex<Double>> f = Fourier::fft(signal.column(k));
        for(UInt i=0; i<f.size(); i++)
          spectrogram(idx+i, 2+k) = std::abs(f.at(i))*std::sqrt(sampling/f.size());
      }
    }, comm);
    Parallel::reduceSum(spectrogram, 0, comm);

    if(Parallel::isMaster(comm))
    {
      logStatus<<"write spectrogram <"<<fileNameOut<<">"<<Log::endl;
      writeFileMatrix(fileNameOut, spectrogram);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
