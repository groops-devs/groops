/***********************************************/
/**
* @file gnssProcessingStepSelectNormalsBlockStructure.h
*
* @brief GNSS processing step: SelectNormalsBlockStructure.
*
* @author Torsten Mayer-Guerr
* @date 2021-09-05
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPROCESSINGSTEPSELECTNORMALSBLOCKSTRUCTURE__
#define __GROOPS_GNSSPROCESSINGSTEPSELECTNORMALSBLOCKSTRUCTURE__

// Latex documentation
#ifdef DOCSTRING_GnssProcessingStep
static const char *docstringGnssProcessingStepSelectNormalsBlockStructure = R"(
\subsection{SelectNormalsBlockStructure}\label{gnssProcessingStepType:selectNormalsBlockStructure}
Select block structure of sparse normal equations for subsequent steps.

This step can be used to define the structure of the different parts of the normal equation system,
which can have a major impact on computing performance and memory consumption depending on the processing setup.

\fig{!hb}{0.4}{gnss_normals_structure}{fig:gnss_normals_structure}{Structure of normal equations in GNSS processing}

The normal equation system is divided into three parts for epoch, interval, and ambiguity parameters.
The epoch part is subdivided further into one subpart per epoch. Each part is divided into blocks and only non-zero
blocks are stored in memory to reduce memory consumption and to prevent unnecessary matrix computations.
\config{defaultBlockSizeEpoch}, \config{defaultBlockSizeInterval}, and \config{defaultBlockSizeAmbiguity} control
the size of the blocks within each part of the normal equations. \config{defaultBlockReceiverCount} can be set to group
a number of receivers into one block within the epoch and interval parts.

If \config{keepEpochNormalsInMemory}=\verb|no| epoch blocks are eliminated after they are set up to reduce the number
of parameters in the normal equation system. \config{defaultBlockCountReduction} controls after how many epoch blocks
an elimination step is performed. For larger processing setups or high sampling rates epoch block elimination is recommended
as the large number of clock parameters require a lot of memory.
)";
#endif

/***********************************************/

#include "config/config.h"
#include "gnss/gnssProcessingStep/gnssProcessingStep.h"

/***** CLASS ***********************************/

/** @brief GNSS processing step: SelectNormalsBlockStructure.
* @ingroup gnssProcessingStepGroup
* @see GnssProcessingStep */
class GnssProcessingStepSelectNormalsBlockStructure : public GnssProcessingStepBase
{
  UInt defaultBlockSizeEpoch, defaultBlockSizeInterval, defaultBlockSizeAmbiguity;
  UInt defaultBlockReceiverCount;
  UInt defaultBlockCountReduction;
  Bool keepEpochNormalsInMemory;
  Bool accumulateEpochObservations;

public:
  GnssProcessingStepSelectNormalsBlockStructure(Config &config);
  void process(GnssProcessingStep::State &state) override;
};

/***********************************************/

inline GnssProcessingStepSelectNormalsBlockStructure::GnssProcessingStepSelectNormalsBlockStructure(Config &config)
{
  try
  {
    readConfig(config, "defaultBlockSizeEpoch",       defaultBlockSizeEpoch,       Config::DEFAULT,  "0",  "block size of epoch parameters, 0: one block");
    readConfig(config, "defaultBlockSizeInterval",    defaultBlockSizeInterval,    Config::DEFAULT,  "64", "block size of interval parameters, 0: one block");
    readConfig(config, "defaultBlockSizeAmbiguity",   defaultBlockSizeAmbiguity,   Config::DEFAULT,  "64", "block size of ambiguity parameters, 0: one block");
    readConfig(config, "defaultBlockReceiverCount",   defaultBlockReceiverCount,   Config::DEFAULT,  "0",  "number of receivers to group into one block for epoch and interval");
    readConfig(config, "defaultBlockCountReduction",  defaultBlockCountReduction,  Config::DEFAULT,  "32", "minimum number of blocks for epoch reduction");
    readConfig(config, "keepEpochNormalsInMemory",    keepEpochNormalsInMemory,    Config::DEFAULT,  "1",  "speeds up processing but uses much more memory");
    readConfig(config, "accumulateEpochObservations", accumulateEpochObservations, Config::DEFAULT,  "0",  "set up all observations per epoch and receiver at once");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssProcessingStepSelectNormalsBlockStructure::process(GnssProcessingStep::State &state)
{
  try
  {
    logStatus<<"=== select block structure =================================="<<Log::endl;
    state.normalEquationInfo.defaultBlockSizeEpoch       = defaultBlockSizeEpoch;
    state.normalEquationInfo.defaultBlockSizeInterval    = defaultBlockSizeInterval;
    state.normalEquationInfo.defaultBlockSizeAmbiguity   = defaultBlockSizeAmbiguity;
    state.normalEquationInfo.defaultBlockReceiverCount   = defaultBlockReceiverCount;
    state.normalEquationInfo.defaultBlockCountReduction  = defaultBlockCountReduction;
    state.normalEquationInfo.keepEpochNormalsInMemory    = keepEpochNormalsInMemory;
    state.normalEquationInfo.accumulateEpochObservations = accumulateEpochObservations;
    logInfo<<"  blockSizeEpoch              = "<<state.normalEquationInfo.defaultBlockSizeEpoch    <<Log::endl;
    logInfo<<"  blockSizeInterval           = "<<state.normalEquationInfo.defaultBlockSizeInterval <<Log::endl;
    logInfo<<"  blockSizeAmbiguity          = "<<state.normalEquationInfo.defaultBlockSizeAmbiguity <<Log::endl;
    logInfo<<"  blockReceiverCount          = "<<state.normalEquationInfo.defaultBlockReceiverCount<<Log::endl;
    logInfo<<"  blockCountReduction         = "<<state.normalEquationInfo.defaultBlockCountReduction<<Log::endl;
    logInfo<<"  keepEpochNormalsInMemory    = "<<(state.normalEquationInfo.keepEpochNormalsInMemory ? "yes" : "no")<<Log::endl;
    logInfo<<"  accumulateEpochObservations = "<<(state.normalEquationInfo.accumulateEpochObservations ? "yes" : "no")<<Log::endl;
    state.changedNormalEquationInfo = TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
