/***********************************************/
/**
* @file instrument2Histogram.cpp
*
* @brief Compute a histogram (arc-wise) from an instrument file
*
* @author Andreas Kvas
* @date 2017-07-09
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program computes the arc-wise histogram from an \file{instrument file}{instrument}.
The output is a \file{matrix}{matrix} with the first column containing the lower bound of each bin.
The other columns contain the histograms for each arc.

\fig{!hb}{0.8}{instrument2Histogram}{fig:instrument2Histogram}{GRACE range-rate residuals of one month (one arc) divided into 50 bins.}
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Compute a histogram from a time series.
* @ingroup programsGroup */
class Instrument2Histogram
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Instrument2Histogram, PARALLEL, "compute a histogram from an instrument file", Instrument, Statistics)
GROOPS_RENAMED_PROGRAM(InstrumentComputeHistogram, Instrument2Histogram, date2time(2020, 7, 7))

/***********************************************/

void Instrument2Histogram::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName         fileNameIn, fileNameOut;
    UInt             selectData;
    ExpressionVariablePtr exprMin, exprMax;
    Bool             relative = TRUE;
    Bool             cumulative = FALSE;
    UInt             binCount = 0;

    readConfig(config, "outputfileMatrix",     fileNameOut, Config::MUSTSET,  "",  "column 1: lower bin bound; columns 2 to N: histogram of each arc");
    readConfig(config, "inputfileInstrument",  fileNameIn,  Config::MUSTSET,  "",  "");
    readConfig(config, "selectDataField",      selectData,  Config::DEFAULT,  "0", "select channel for histogram computation");
    readConfig(config, "binCount",             binCount,    Config::OPTIONAL, "",  "(default: Freedman-Diaconis' choice, maximum of all channels)");
    readConfig(config, "lowerBound",           exprMin,     Config::DEFAULT,  "dataMin", "lower bound for bins (default: global minimum, data values outside are ignored)");
    readConfig(config, "upperBound",           exprMax,     Config::DEFAULT,  "dataMax", "upper bound for bins (default: global maximum, data values outside are ignored)");
    readConfig(config, "relative",             relative,    Config::DEFAULT,  "1", "output relative frequencies");
    readConfig(config, "cumulative",           cumulative,  Config::DEFAULT,  "0", "accumulate frequencies");
    if(isCreateSchema(config)) return;

    logStatus<<"reading instrument file <"<<fileNameIn<<">"<<Log::endl;
    InstrumentFile instrumentFile(fileNameIn);
    const UInt arcCount = instrumentFile.arcCount();

    // collect all data
    // ----------------
    std::vector<std::vector<Double>> arcWiseData(arcCount);
    Parallel::forEach(arcWiseData, [&](UInt arcNo)
    {
      Matrix data = instrumentFile.readArc(arcNo).matrix();
      std::vector<Double> values(data.rows());
      for(UInt k=0; k<data.rows(); k++)
        values.at(k) = data(k, selectData+1);
      return values;
    }, comm);
    Parallel::broadCast(arcWiseData, 0, comm);

    // determine bins
    // --------------
    std::vector<Double> bins;
    if(Parallel::isMaster(comm))
    {
      std::vector<Double> globalData;
      for(UInt arcNo=0; arcNo<arcWiseData.size(); arcNo++)
        globalData.insert(globalData.end(), arcWiseData[arcNo].begin(), arcWiseData[arcNo].end());
      std::sort(globalData.begin(), globalData.end());

      auto varList = config.getVarList();
      addVariable("dataMin", globalData.front(), varList);
      addVariable("dataMax", globalData.back(),  varList);
      const Double lowerBound = exprMin->evaluate(varList);
      const Double upperBound = exprMax->evaluate(varList);

      const UInt orginalSize = globalData.size();
      globalData.erase(std::remove_if(globalData.begin(), globalData.end(), [lowerBound, upperBound](Double x) {return (x<lowerBound) || (x>upperBound);}), globalData.end());
      logInfo<<"Discarded "<<orginalSize-globalData.size()<<" elements."<<Log::endl;

      // compute number of bins based on Freedman-Diaconis' choice
      if(binCount == 0)
      {
        const UInt   count       = (globalData.size()+1)/2;
        const Double q1          = globalData[count / 2];
        const Double q3          = globalData[3 * count / 2];
        const Double binSize     = 2.0 * (q3 - q1) / std::pow(globalData.size(), 1./3.);
        const UInt   binCountMax = 100;
        binCount = std::min(binCountMax, static_cast<UInt>(std::ceil((upperBound-lowerBound)/binSize)));
        if(binCount == binCountMax)
          logWarning << "Bin count set to a maximum of <" << binCountMax << ">" << Log::endl;
      }
      logInfo<<"Sort data into "<<binCount<<" bins in the range of ["<<lowerBound<<", "<<upperBound<<"]"<<Log::endl;

      bins.resize(1, lowerBound); // intervals
      for(UInt k=0; k<binCount; k++)
        bins.push_back(bins.back() + (upperBound-lowerBound)/binCount);
      bins.back() = upperBound; // make sure upper bound is correct
    }
    Parallel::broadCast(bins, 0, comm);

    // compute histogram
    // -----------------
    logStatus<<"compute histogram"<<Log::endl;
    Matrix histogram(bins.size()-1, arcCount+1); // first column: lower bin bound

    Parallel::forEach(arcCount, [&](UInt arcNo)
    {
      std::vector<Double> data = arcWiseData.at(arcNo);
      for(UInt k = 0; k<bins.size()-2; k++)
        histogram(k, arcNo+1) = std::count_if(data.begin(), data.end(), [&, bins, k](Double v){ return (v>=bins[k] && v<bins[k+1]); });

      // last bin includes upper bound
      histogram(bins.size()-2, arcNo+1) = std::count_if(data.begin(), data.end(), [&, bins](Double v){ return (v>=bins[bins.size()-2] && v<=bins.back()); });

      if(relative)
      {
        UInt count = std::count_if(data.begin(), data.end(), [&, bins](Double v){ return (v>=bins.front() && v<=bins.back()); });
        histogram.column(arcNo+1) *= 1./count;
      }
      if(cumulative)
      {
        for(UInt k = 1; k<bins.size()-1; k++)
          histogram(k, arcNo+1) += histogram(k-1, arcNo+1);
      }
    }, comm);
    Parallel::reduceSum(histogram, 0, comm);

    logStatus<<"write histogram to <"<<fileNameOut<<">"<<Log::endl;
    if(Parallel::isMaster(comm))
    {
      for(UInt k = 0; k<bins.size()-1; k++)
        histogram(k, 0) = bins.at(k);
      writeFileMatrix(fileNameOut, histogram);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
