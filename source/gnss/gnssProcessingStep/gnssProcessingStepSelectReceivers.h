/***********************************************/
/**
* @file gnssProcessingStepSelectReceivers.h
*
* @brief GNSS processing step: SelectReceivers.
*
* @author Torsten Mayer-Guerr
* @date 2021-09-05
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPROCESSINGSTEPSELECTRECEIVERS__
#define __GROOPS_GNSSPROCESSINGSTEPSELECTRECEIVERS__

// Latex documentation
#ifdef DOCSTRING_GnssProcessingStep
static const char *docstringGnssProcessingStepSelectReceivers = R"(
\subsection{SelectReceivers}\label{gnssProcessingStepType:selectReceivers}
This step can be used to process only a subset of stations in subsequent processing steps.
The most common use is to start the processing with a well-distributed network of core stations as seen in
\reference{GNSS satellite orbit determination and network analysis}{cookbook.gnssNetwork:processing}.
To later process all other stations individually, use the processing step
\configClass{processingStep:forEachReceiverSeparately}{gnssProcessingStepType:forEachReceiverSeparately}
and select all stations excluding the core stations in that step.
)";
#endif

/***********************************************/

#include "config/config.h"
#include "classes/platformSelector/platformSelector.h"
#include "gnss/gnssProcessingStep/gnssProcessingStep.h"

/***** CLASS ***********************************/

/** @brief GNSS processing step: SelectReceivers.
* @ingroup gnssProcessingStepGroup
* @see GnssProcessingStep */
class GnssProcessingStepSelectReceivers : public GnssProcessingStepBase
{
  PlatformSelectorPtr selectReceivers;

public:
  GnssProcessingStepSelectReceivers(Config &config);
  void process(GnssProcessingStep::State &state) override;
  Bool expectInitializedParameters() const override {return FALSE;}
};

/***********************************************/

inline GnssProcessingStepSelectReceivers::GnssProcessingStepSelectReceivers(Config &config)
{
  try
  {
    readConfig(config, "selectReceivers", selectReceivers, Config::MUSTSET,  "", "");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void GnssProcessingStepSelectReceivers::process(GnssProcessingStep::State &state)
{
  try
  {
    logStatus<<"=== select receivers ========================================"<<Log::endl;
    if(state.normalEquationInfo.isEachReceiverSeparately)
    {
      logWarning<<"SelectReceivers is not allowed in single receiver loop"<<Log::endl;
      return;
    }
    state.normalEquationInfo.estimateReceiver = state.gnss->selectReceivers(selectReceivers);
    for(auto recv : state.gnss->receivers)
      if(!recv->useable())
        state.normalEquationInfo.estimateReceiver.at(recv->idRecv()) = FALSE;
    logInfo<<"  "<<std::count(state.normalEquationInfo.estimateReceiver.begin(), state.normalEquationInfo.estimateReceiver.end(), TRUE)<<" receivers selected"<<Log::endl;
    state.changedNormalEquationInfo = TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
