/***********************************************/
/**
* @file gnssProcessingStepSelectEpochs.h
*
* @brief GNSS processing step: SelectEpochs.
*
* @author Torsten Mayer-Guerr
* @date 2021-09-05
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPROCESSINGSTEPSELECTEPOCHS__
#define __GROOPS_GNSSPROCESSINGSTEPSELECTEPOCHS__

// Latex documentation
#ifdef DOCSTRING_GnssProcessingStep
static const char *docstringGnssProcessingStepSelectEpochs = R"(
\subsection{SelectEpochs}\label{gnssProcessingStepType:selectEpochs}
Select epochs for subsequent steps. This step can be used to reduce the processing sampling
while keeping the original observation sampling for all preprocessing steps (e.g. outlier and cycle slip detection).
Another example is to process at a 5-minute sampling by setting \config{nthEpoch}=\verb|10| and then
at the end to densify only the clock parameters to the full 30-second observation sampling by
setting \config{nthEpoch}=\verb|1| while keeping all other parameters fixed
with \configClass{processingStep:selectParametrizations}{gnssProcessingStepType:selectParametrizations}.
)";
#endif

/***********************************************/

#include "config/config.h"
#include "gnss/gnssProcessingStep/gnssProcessingStep.h"

/***** CLASS ***********************************/

/** @brief GNSS processing step: SelectEpochs.
* @ingroup gnssProcessingStepGroup
* @see GnssProcessingStep */
class GnssProcessingStepSelectEpochs : public GnssProcessingStepBase
{
  UInt nthEpoch;

public:
  GnssProcessingStepSelectEpochs(Config &config);
  void process(GnssProcessingStep::State &state) override;
  Bool expectInitializedParameters() const override {return FALSE;}
};

/***********************************************/

inline GnssProcessingStepSelectEpochs::GnssProcessingStepSelectEpochs(Config &config)
{
  try
  {
    readConfig(config, "nthEpoch", nthEpoch, Config::MUSTSET, "1", "use only every nth epoch in all subsequent processing steps");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssProcessingStepSelectEpochs::process(GnssProcessingStep::State &state)
{
  try
  {
    logStatus<<"=== select epochs ==========================================="<<Log::endl;
    state.normalEquationInfo.idEpochs.clear();
    for(UInt idEpoch=0; idEpoch<state.gnss->times.size(); idEpoch+=nthEpoch)
      state.normalEquationInfo.idEpochs.push_back(idEpoch);
    logInfo<<"  "<<state.normalEquationInfo.idEpochs.size()<<" of "<<state.gnss->times.size()<<" epochs selected"<<Log::endl;
    state.changedNormalEquationInfo = TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
