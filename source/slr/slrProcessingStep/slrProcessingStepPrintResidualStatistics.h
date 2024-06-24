/***********************************************/
/**
* @file slrProcessingStepPrintResidualStatistics.h
*
* @brief SLR processing step: PrintResidualStatistics.
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#ifndef __GROOPS_SLRPROCESSINGSTEPPRINTRESIDUALSTATISTICS__
#define __GROOPS_SLRPROCESSINGSTEPPRINTRESIDUALSTATISTICS__

// Latex documentation
#ifdef DOCSTRING_SlrProcessingStep
static const char *docstringSlrProcessingStepPrintResidualStatistics = R"(
\subsection{PrintResidualStatistics}\label{slrProcessingStepType:printResidualStatistics}
Print residual statistics.
\begin{verbatim}
  station   sigma redundancy obsCount outlier
  ----------------------------------------------
  1874       0.52   0.86      22      1 (4.55 %)
  1889       1.20   0.98     186      5 (2.69 %)
  1890       0.63   0.77      14      1 (7.14 %)
  1891       0.49   0.50       6      0 (0.00 %)
  7237       1.08   0.95     236     14 (5.93 %)
  7394       0.36   0.88      26      0 (0.00 %)
  7811       0.38   0.41       5      0 (0.00 %)
  7819       1.21   0.94     120      1 (0.83 %)
  7821       0.69   0.95     202      3 (1.49 %)
  7827       0.40   0.85      29      1 (3.45 %)
  7839       0.52   0.93     143     10 (6.99 %)
  7840       0.15   0.80      16      0 (0.00 %)
  7841       0.26   0.90      56      1 (1.79 %)
  7941       0.55   0.92     277      5 (1.81 %)
  8834       0.66   0.88     101      1 (0.99 %)
  ----------------------------------------------
  satellite sigma redundancy obsCount outlier
  ----------------------------------------------
  lageos1    1.04   0.94     722     24 (3.32 %)
  lageos2    0.91   0.95     590     11 (1.86 %)
  etalon1    1.19   0.78      57      2 (3.51 %)
  etalon2    1.10   0.81      70      6 (8.57 %)
  ----------------------------------------------
\end{verbatim}
)";
#endif

/***********************************************/

#include "config/config.h"
#include "slr/slrProcessingStep/slrProcessingStep.h"

/***** CLASS ***********************************/

/** @brief SLR processing step: PrintResidualStatistics.
* @ingroup slrProcessingStepGroup
* @see SlrProcessingStep */
class SlrProcessingStepPrintResidualStatistics : public SlrProcessingStepBase
{
public:
  SlrProcessingStepPrintResidualStatistics(Config &/*config*/) {}
  void process(SlrProcessingStep::State &state) override;
};

/***********************************************/

inline void SlrProcessingStepPrintResidualStatistics::process(SlrProcessingStep::State &state)
{
  try
  {
    logStatus<<"=== print residual statistics  =============================="<<Log::endl;
    constexpr Double huber      = 2.5;
    constexpr Double huberPower = 1.5;

    logInfo<<"  station   sigma redundancy obsCount outlier"<<Log::endl;
    logInfo<<"  ----------------------------------------------"<<Log::endl;
    for(UInt idStat=0; idStat<state.slr->stations.size(); idStat++)
      if(state.normalEquationInfo.estimateStation.at(idStat))
      {
        Double  ePe = 0, redundancy = 0;
        UInt    obsCount = 0, outlierCount = 0;
        state.residualsStatistics(idStat, NULLINDEX/*idSat*/, ePe, redundancy, obsCount, outlierCount);
        if(obsCount)
        {
          std::string name = state.slr->stations.at(idStat)->name();
          name.resize(10, ' ');
          logInfo<<"  "<<name.substr(0, 10)
                 <<state.sigmaFactor(idStat) * Vce::standardDeviation(ePe, redundancy, huber, huberPower)%"%5.2f"s<<(redundancy/obsCount)%"%7.2f"s
                 <<obsCount%"%8i"s<<outlierCount%"%7i ("s<<(100.*outlierCount/obsCount)%"%4.2f %%)"s<<Log::endl;
        }
      }

    logInfo<<"  ----------------------------------------------"<<Log::endl;
    logInfo<<"  satellite sigma redundancy obsCount outlier"<<Log::endl;
    logInfo<<"  ----------------------------------------------"<<Log::endl;
    for(UInt idSat=0; idSat<state.slr->satellites.size(); idSat++)
      if(state.normalEquationInfo.estimateSatellite.at(idSat))
      {
        Double  ePe = 0, redundancy = 0;
        UInt    obsCount = 0, outlierCount = 0;
        state.residualsStatistics(NULLINDEX/*idStat*/, idSat, ePe, redundancy, obsCount, outlierCount);
        if(obsCount)
        {
          std::string name = state.slr->satellites.at(idSat)->name();
          name.resize(10, ' ');
          logInfo<<"  "<<name.substr(0, 10)
                 <<Vce::standardDeviation(ePe, redundancy, huber, huberPower)%"%5.2f"s<<(redundancy/obsCount)%"%7.2f"s
                 <<obsCount%"%8i"s<<outlierCount%"%7i ("s<<(100.*outlierCount/obsCount)%"%4.2f %%)"s<<Log::endl;
        }
      }
    logInfo<<"  ----------------------------------------------"<<Log::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
