/***********************************************/
/**
* @file gnssProcessingStepPrintResidualStatistics.h
*
* @brief GNSS processing step: PrintResidualStatistics.
*
* @author Torsten Mayer-Guerr
* @date 2021-09-05
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPROCESSINGSTEPPRINTRESIDUALSTATISTICS__
#define __GROOPS_GNSSPROCESSINGSTEPPRINTRESIDUALSTATISTICS__

// Latex documentation
#ifdef DOCSTRING_GnssProcessingStep
static const char *docstringGnssProcessingStepPrintResidualStatistics = R"(
\subsection{PrintResidualStatistics}\label{gnssProcessingStepType:printResidualStatistics}
Print residual statistics.
\begin{verbatim}
  areq: C1CG**: factor =  0.64, sigma0 = 1.00, count =  2748, outliers =    48 (1.75 \%)
  areq: C1WG**: factor =  0.50, sigma0 = 1.00, count =  2748, outliers =    43 (1.56 \%)
  areq: C2WG**: factor =  0.50, sigma0 = 1.00, count =  2748, outliers =    59 (2.15 \%)
  areq: C5XG**: factor =  0.46, sigma0 = 1.00, count =  1279, outliers =    23 (1.80 \%)
  areq: L1CG**: factor =  0.86, sigma0 = 0.96, count =  2748, outliers =    40 (1.46 \%)
  areq: L1WG**: factor =  0.86, sigma0 = 1.02, count =  2748, outliers =    40 (1.46 \%)
  areq: L2WG**: factor =  0.86, sigma0 = 0.96, count =  2748, outliers =    40 (1.46 \%)
  areq: L5XG**: factor =  0.86, sigma0 = 1.30, count =  1279, outliers =    14 (1.09 \%)
  areq: C1PR**: factor =  0.48, sigma0 = 1.00, count =  1713, outliers =    53 (3.09 \%)
  areq: C2PR**: factor =  0.55, sigma0 = 1.00, count =  1713, outliers =    51 (2.98 \%)
  areq: L1PR**: factor =  0.85, sigma0 = 1.09, count =  1713, outliers =    29 (1.69 \%)
  areq: L2PR**: factor =  0.85, sigma0 = 0.88, count =  1713, outliers =    29 (1.69 \%)
  areq: C1XE**: factor =  0.44, sigma0 = 1.00, count =  1264, outliers =    21 (1.66 \%)
  areq: C5XE**: factor =  0.33, sigma0 = 1.00, count =  1264, outliers =    27 (2.14 \%)
  areq: C7XE**: factor =  0.28, sigma0 = 1.00, count =  1264, outliers =    41 (3.24 \%)
  areq: L1XE**: factor =  0.82, sigma0 = 1.14, count =  1264, outliers =    15 (1.19 \%)
  areq: L5XE**: factor =  0.82, sigma0 = 0.84, count =  1264, outliers =    15 (1.19 \%)
  areq: L7XE**: factor =  0.82, sigma0 = 0.94, count =  1264, outliers =    15 (1.19 \%)
  badg: C1CG**: factor =  1.25, sigma0 = 1.00, count =  2564, outliers =    47 (1.83 \%)
  ...
\end{verbatim}
)";
#endif

/***********************************************/

#include "config/config.h"
#include "gnss/gnssProcessingStep/gnssProcessingStep.h"

/***** CLASS ***********************************/

/** @brief GNSS processing step: PrintResidualStatistics.
* @ingroup gnssProcessingStepGroup
* @see GnssProcessingStep */
class GnssProcessingStepPrintResidualStatistics : public GnssProcessingStepBase
{
public:
  GnssProcessingStepPrintResidualStatistics(Config &/*config*/) {}
  void process(GnssProcessingStep::State &state) override;
};

/***********************************************/

inline void GnssProcessingStepPrintResidualStatistics::process(GnssProcessingStep::State &state)
{
  try
  {
    logStatus<<"=== print residual statistics  =============================="<<Log::endl;
    for(UInt idRecv=0; idRecv<state.gnss->receivers.size(); idRecv++)
      if(state.normalEquationInfo.estimateReceiver.at(idRecv))
      {
        std::vector<GnssType> types = state.gnss->types(~(GnssType::PRN + GnssType::FREQ_NO));
        std::vector<Double>   ePe(types.size(), 0), redundancy(types.size(), 0);
        std::vector<UInt>     obsCount(types.size(), 0), outlierCount(types.size(), 0);
        state.residualsStatistics(idRecv, NULLINDEX/*idTrans*/, types, ePe, redundancy, obsCount, outlierCount);
        Vector factors(types.size());
        for(UInt i=0; i<types.size(); i++)
        {
          const UInt idx = GnssType::index(state.stations.at(idRecv).sigmaTypes, types.at(i));
          if(idx != NULLINDEX)
            factors(i) = state.stations.at(idRecv).sigmaFactors.at(idx);
        }
        Parallel::reduceSum(factors, 0, state.normalEquationInfo.comm);
        if(Parallel::isMaster(state.normalEquationInfo.comm))
          for(UInt i=0; i<types.size(); i++)
            if(obsCount.at(i))
            {
              logInfo<<"  "<<state.gnss->receivers.at(idRecv)->name()<<": "<<types.at(i).str()
                    <<": factor = "    <<factors(i)%"%5.2f"s
                    <<", sigma0 = "    <<Vce::standardDeviation(ePe.at(i), redundancy.at(i), 2.5/*huber*/, 1.5/*huberPower*/)%"%4.2f"s
                    <<", count = "     <<obsCount.at(i)%"%5i"s
                    <<", outliers = "  <<outlierCount.at(i)%"%5i"s<<" ("<<(100.*outlierCount.at(i)/obsCount.at(i))%"%4.2f"s<<" %)"<<Log::endl;
            }
      } // for(idRecv)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
