/***********************************************/
/**
* @file instrumentWaveletDecomposition.cpp
*
* @brief Wavelet decomposition of an instrument file.
*
* @author Saniya Behzadpour
* @date 2018-07-17
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program performs a multilevel one-dimensional wavelet analysis on one \config{selectDataField}
data column of \configFile{inputfileInstrument}{instrument}.
The \configFile{outputfileInstrument}{instrument} contains the decomposed levels in time domain ${a_J,d_J,...,d_1}$
)";


/***********************************************/
#include "base/wavelets.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"
#include "parallel/parallel.h"
#include "programs/program.h"

/***** CLASS ***********************************/

/** @brief Wavelet decomposition of an instrument file.
* @ingroup programsGroup */
class InstrumentWaveletDecomposition
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(InstrumentWaveletDecomposition, SINGLEPROCESS, "Wavelet decomposition of an instrument file.", Instrument)

/***********************************************/

void InstrumentWaveletDecomposition::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameOut;
    FileName fileNameIn, fileNameInWavelet;
    UInt     selectData;
    UInt     level;

    readConfig(config, "outputfileInstrument", fileNameOut,       Config::MUSTSET, "", "MISCVALUES, decomposed levels in time domain a_J,d_J,...,d_1");
    readConfig(config, "inputfileInstrument",  fileNameIn,        Config::MUSTSET, "", "");
    readConfig(config, "selectDataField",      selectData,        Config::DEFAULT, "0", "select a data column for decomposition");
    readConfig(config, "inputfileWavelet",     fileNameInWavelet, Config::MUSTSET, "{groopsDataDir}/wavelets/", "wavelet coefficients");
    readConfig(config, "level",                level,             Config::MUSTSET, "4", "level of decomposition");
    if(isCreateSchema(config)) return;

    logStatus<<"read instrument file "<<"<"<<fileNameIn<<">"<<Log::endl;
    const Arc arc = InstrumentFile::read(fileNameIn);

    logStatus<<"read wavelet filters "<<"<"<<fileNameInWavelet<<">"<<Log::endl;
    Vector wl;
    readFileMatrix(fileNameInWavelet, wl);

    // perform n-Level DWT
    std::vector<UInt>   length, flag;
    std::vector<Double> dwtCoeff;
    std::vector<Double> signal = Vector(arc.matrix().column(1+selectData));
    Wavelets::discreteWaveletTransform(signal, wl, level, dwtCoeff, flag, length);

    // detail coefficients in each column
    UInt idx = 0;
    std::vector<std::vector<Double>> dwt(level+1, std::vector<Double>(dwtCoeff.size(), 0.));
    for(UInt j=0; j<level+1;j++)
      for(UInt i=0; i<length.at(j); i++)
      {
        dwt.at(j).at(idx) = dwtCoeff.at(idx);
        idx++;
      }

    // apply inverse DWT, from wavelet domain to time domain
    std::vector<std::vector<Double>> idwt(level+1);
    for (UInt j=0; j<level+1; j++)
      Wavelets::inverseDiscreteWaveletTransform(dwt.at(j), wl, flag, idwt.at(j), length);

    MiscValuesArc arcOut;
    for(UInt i=0; i<arc.size(); i++)
    {
      MiscValuesEpoch epoch(level+1);
      epoch.time = arc.at(i).time;
      for(UInt j=0; j<level+1; j++)
        epoch.values(j) = idwt.at(j).at(i);
      arcOut.push_back(epoch);
    }

    logStatus<<"write decomposed signals"<<"<"<<fileNameOut<<">"<<Log::endl;
    InstrumentFile::write(fileNameOut, arcOut);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/*************************************************/
