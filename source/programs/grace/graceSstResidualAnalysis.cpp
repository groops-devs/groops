/***********************************************/
/**
* @file graceSstResidualAnalysis.cpp
*
* @brief Multiresolution analysis of GRACE SST residuals.
*
* @author Saniya Behzadpour
* @date 2018-07-15
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program applies the Multi-Resolution Analysis (MRA) using
Discrete Wavelet Transform (DWT) to the monthly GRACE SST post-fit residuals.
First, the residuals are transfered into wavelet domain by applying an 8 level
Daubechies wavelet transform (default).
In the next step, detail coefficients are merged into three major groups
due to their approximate frequency subbands:
\begin{itemize}
\item Low scale details, corresponding to the frequency band above 10 mHz;
\item Intermediate scale details, corresponding to the approximate frequency
      range above 3 mHz up to 10 mHz;
\item High scale details, corresponding to the approximate frequency range
above 0.5 mHz up to 10 mHz.
\end{itemize}
In the last step, each group is reconstructed back into time domain.
)";

/***********************************************/

#include "base/wavelets.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"
#include "parallel/parallel.h"
#include "programs/program.h"

/***** CLASS ***********************************/

/** @brief Multiresolution analysis of GRACE SST residuals.
* @ingroup programsGroup */
class GraceSstResidualAnalysis
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GraceSstResidualAnalysis, SINGLEPROCESS, "Multiresolution analysis of GRACE SST residuals.", Grace)

/***********************************************/

void GraceSstResidualAnalysis::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameOutHigh, fileNameOutMid, fileNameOutLow;
    FileName fileNameInResiduals, fileNameInWavelet;
    UInt     level = 8;

    readConfig(config, "outputfileInstrumentHighScale", fileNameOutHigh,     Config::MUSTSET, "", "High scale details");
    readConfig(config, "outputfileInstrumentMidScale",  fileNameOutMid,      Config::MUSTSET, "", "Intermediate scale details");
    readConfig(config, "outputfileInstrumentLowScale",  fileNameOutLow,      Config::MUSTSET, "", "Low scale details");
    readConfig(config, "inputfileInstrument",           fileNameInResiduals, Config::MUSTSET, "", "GRACE SST Residuals");
    readConfig(config, "inputfileWavelet",              fileNameInWavelet,   Config::MUSTSET, "{groopsDataDir}/wavelets/db20.txt", "wavelet coefficients");
    if(isCreateSchema(config)) return;

    // =======================

    logStatus<<"read GRACE SST residuals "<<"<"<<fileNameInResiduals<<">"<<Log::endl;
    SatelliteTrackingArc sstArc = InstrumentFile::read(fileNameInResiduals);
    UInt posCount = sstArc.size();
    std::vector<Double> signal(posCount);
    for(UInt i=0; i<posCount; i++)
      signal.at(i) = sstArc.at(i).rangeRate;

    logStatus<<"read wavelet filters "<<"<"<<fileNameInWavelet<<">"<<Log::endl;
    Vector wl;
    readFileMatrix(fileNameInWavelet, wl);

    std::vector<UInt>   length;
    std::vector<UInt>   flag, flag2;
    std::vector<Double> dwtCoeff,dwt_a08,idwt_a08,idwt_out;
    std::vector<Double> dwt_d03, dwt_d02, dwt_d01;
    std::vector<Double> idwt_d03,idwt_d02,idwt_d01;

    //perform 8-Level DWT
    Wavelets::discreteWaveletTransform(signal, wl, level, dwtCoeff,flag, length);

    //compute cumulative flag vector
    flag2.resize(length.size());
    for (UInt i=0; i<length.size();i++)
      flag2.at(i) = length.at(i);
    for (UInt i=1; i<flag2.size();i++)
      flag2.at(i) += flag2.at(i-1);

    //merge detail coefficients into three major groups
    dwt_a08.resize(dwtCoeff.size()); dwt_d03.resize(dwtCoeff.size());
    dwt_d02.resize(dwtCoeff.size()); dwt_d01.resize(dwtCoeff.size());
    for (UInt i = 0; i < dwtCoeff.size();i++)
      if (i < flag2.at(0))
        dwt_a08.at(i) = dwtCoeff.at(i);
      else if (i < flag2.at(3))
        dwt_d03.at(i) = dwtCoeff.at(i);
      else if (i < flag2.at(5))
        dwt_d02.at(i) = dwtCoeff.at(i);
      else
        dwt_d01.at(i) = dwtCoeff.at(i);

    //apply inverse DWT, from wavelet domain to time domain
    Wavelets::inverseDiscreteWaveletTransform(dwt_d03, wl, flag, idwt_d03, length);
    Wavelets::inverseDiscreteWaveletTransform(dwt_d02, wl, flag, idwt_d02, length);
    Wavelets::inverseDiscreteWaveletTransform(dwt_d01, wl, flag, idwt_d01, length);

    SatelliteTrackingArc d03, d02, d01;
    for(UInt i=0; i<posCount; i++)
    {
      SatelliteTrackingEpoch epoch;
      epoch.time = sstArc.at(i).time;
      epoch.range = epoch.rangeRate = epoch.rangeAcceleration = idwt_d03.at(i);
      d03.push_back(epoch);
    }

    for(UInt i=0; i<posCount; i++)
    {
      SatelliteTrackingEpoch epoch;
      epoch.time = sstArc.at(i).time;
      epoch.range = epoch.rangeRate = epoch.rangeAcceleration = idwt_d02.at(i);
      d02.push_back(epoch);
    }

    for(UInt i=0; i<posCount; i++)
    {
      SatelliteTrackingEpoch epoch;
      epoch.time = sstArc.at(i).time;
      epoch.range = epoch.rangeRate = epoch.rangeAcceleration = idwt_d01.at(i);
      d01.push_back(epoch);
    }

    logStatus<<"write decomposed signals"<<Log::endl;
    logStatus<<"<"<<fileNameOutHigh<<">"<<Log::endl;
    InstrumentFile::write(fileNameOutHigh, d03);
    logStatus<<"<"<<fileNameOutMid<<">"<<Log::endl;
    InstrumentFile::write(fileNameOutMid,  d02);
    logStatus<<"<"<<fileNameOutLow<<">"<<Log::endl;
    InstrumentFile::write(fileNameOutLow,  d01);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}
/*************************************************/
