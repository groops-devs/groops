/***********************************************/
/**
* @file instrument2Scaleogram.cpp
*
* @brief Compute scalogram and detail levels of an instrument file.
*
* @author Andreas Kvas
* @author Saniya Behzadpour
* @date 2018-05-01
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program computes the wavelet transform of a time series up to a \config{maxLevel}.
The scalogram is written to a matrix which can be plotted by using a gridded layer in \program{PlotGraph}.
Individual detail levels can be written to matrix files by setting \configFile{outputfileLevels}{matrix}.
The data column to be decomposed must be set by \config{selectDataField}.

The wavelet transform is implemented as a filter bank, so care should be taken when the input contains data gaps.
Low/highpass wavelet filters are applied in forward and backward direction, input is padded symmetric.
See \configClass{digitalFilter}{digitalFilterType} for details.

\fig{!hb}{0.8}{instrument2Scaleogram}{fig:Instrument2Scaleogram}{GRACE range-rate residuals of one month.}
)";

/***********************************************/

#include "programs/program.h"
#include "base/wavelets.h"
#include "base/polynomial.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"


/***** CLASS ***********************************/

/** @brief Compute scalogram and detail levels of an instrument file.
* @ingroup programsGroup */
class Instrument2Scaleogram
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Instrument2Scaleogram, SINGLEPROCESS, "Compute scalogram and detail levels of an instrument file.", Instrument, Statistics)
GROOPS_RENAMED_PROGRAM(InstrumentComputeScaleogram, Instrument2Scaleogram, date2time(2020, 7, 7))

/***********************************************/

void Instrument2Scaleogram::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameIn, fileNameOutGrid, fileNameLevels;
    FileName fileNameWavelet;
    UInt     dataField;
    UInt     maxLevel = MAX_UINT;

    readConfig(config, "outputfileScaleogram", fileNameOutGrid, Config::MUSTSET,  "",  "matrix columns: mjd, level, value");
    readConfig(config, "outputfileLevels",     fileNameLevels,  Config::OPTIONAL, "",  "use loopLevel as variable");
    readConfig(config, "inputfileInstrument",  fileNameIn,      Config::MUSTSET,  "",  "");
    readConfig(config, "inputfileWavelet",     fileNameWavelet, Config::MUSTSET,  "{groopsDataDir}/wavelets/", "wavelet coefficients");
    readConfig(config, "selectDataField",      dataField,       Config::DEFAULT,  "0", "data column to transform");
    readConfig(config, "maxLevel",             maxLevel,        Config::OPTIONAL, "",  "maximum level of decomposition (default: full)");
    if(isCreateSchema(config)) return;

    VariableList fileNameVariableList;
    addVariable("loopLevel", fileNameVariableList);

    logStatus<<"read instrument file <"<<fileNameIn<<">"<<Log::endl;
    Arc arc = InstrumentFile::read(fileNameIn);
    std::vector<Time> times = arc.times();
    Matrix data = arc.matrix();

    Vector wl;
    readFileMatrix(fileNameWavelet, wl);
    wl *= 1./std::sqrt(2); // normalized wavelet

    logStatus<<"compute wavelet transform"<<Log::endl;
    std::vector<Matrix> levels = Wavelets::waveletTransform(data.column(dataField+1), wl, maxLevel);
    logInfo<<"  levels = "<<levels.size()<<Log::endl;

    Matrix scaleogramGrid(times.size()*levels.size(), 3);
    Polynomial poly(0); // interpolation
    for(UInt k=0; k<levels.size(); k++)
    {
      std::vector<Time> interpTimes(levels.at(k).rows());
      for(UInt i=0; i<interpTimes.size(); i++)
        interpTimes.at(i) = times.front() + 1.*i/(interpTimes.size()-1) * (times.back()-times.front());

      if(!fileNameLevels.empty())
      {
        Matrix output(levels.at(k).rows(), 2);
        for(UInt i=0; i<interpTimes.size(); i++)
        {
          output(i, 0) = interpTimes.at(i).mjd();
          output(i, 1) = levels.at(k)(i,0);
        }
        fileNameVariableList["loopLevel"]->setValue(k);
        writeFileMatrix(fileNameLevels(fileNameVariableList), output);
      }

      copy(data.column(0),                                     scaleogramGrid.slice(k*times.size(), 0, times.size(), 1));
      copy(Vector(data.rows(), k),                             scaleogramGrid.slice(k*times.size(), 1, times.size(), 1));
      copy(poly.interpolate(times, interpTimes, levels.at(k)), scaleogramGrid.slice(k*times.size(), 2, times.size(), 1));
    }

    logStatus<<"write scaleogram to <"<<fileNameOutGrid<<">"<<Log::endl;
    writeFileMatrix(fileNameOutGrid, scaleogramGrid);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
