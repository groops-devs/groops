/***********************************************/
/**
* @file gnssProcessingStepForEachReceiverSeparately.h
*
* @brief GNSS processing step: ForEachReceiverSeparately.
*
* @author Torsten Mayer-Guerr
* @date 2021-09-05
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPROCESSINGSTEPFOREACHRECEIVERSEPARATELY__
#define __GROOPS_GNSSPROCESSINGSTEPFOREACHRECEIVERSEPARATELY__

// Latex documentation
#ifdef DOCSTRING_GnssProcessingStep
static const char *docstringGnssProcessingStepForEachReceiverSeparately = R"(
\subsection{ForEachReceiverSeparately}\label{gnssProcessingStepType:forEachReceiverSeparately}
Perform these processing steps for each \configClass{selectReceivers}{gnssTransceiverSelectorType} separately.
All non-receiver related parameters parameters are disabled in these processing steps (see .

This step can be used for individual precise point positioning (PPP) of all stations.
During \reference{GNSS satellite orbit determination and network analysis}{cookbook.gnssNetwork:processing} this step is used after the
initial processing of the core network to process all other stations individually. In that case provide the same station list as
\configFile{inputfileExcludeStationList}{stringList} in this step that was used as \configFile{inputfileStationList}{stringList} in the
\configClass{selectReceivers}{gnssProcessingStepType:selectReceivers} step where the core network was selected.
)";
#endif

/***********************************************/

#include "config/config.h"
#include "gnss/gnssTransceiverSelector/gnssTransceiverSelector.h"
#include "gnss/gnssProcessingStep/gnssProcessingStep.h"

/***** CLASS ***********************************/

/** @brief GNSS processing step: ForEachReceiverSeparately.
* @ingroup gnssProcessingStepGroup
* @see GnssProcessingStep */
class GnssProcessingStepForEachReceiverSeparately : public GnssProcessingStepBase
{
  GnssTransceiverSelectorPtr selectReceivers;
  std::string                variableReceiver;
  Config                     configProcessingSteps;

public:
  GnssProcessingStepForEachReceiverSeparately(Config &config);
  void process(GnssProcessingStep::State &state) override;
  Bool expectInitializedParameters() const override {return FALSE;}
};

/***********************************************/

inline GnssProcessingStepForEachReceiverSeparately::GnssProcessingStepForEachReceiverSeparately(Config &config)
{
  try
  {
    GnssProcessingStepPtr processingSteps;

    readConfig(config, "selectReceivers",  selectReceivers,  Config::MUSTSET,  "",        "");
    readConfig(config, "variableReceiver", variableReceiver, Config::OPTIONAL, "station", "variable is set for each receiver");
    readConfigLater(config, "processingStep", processingSteps, configProcessingSteps, Config::MUSTSET, "", "steps are processed consecutively");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void GnssProcessingStepForEachReceiverSeparately::process(GnssProcessingStep::State &state)
{
  try
  {
    Parallel::barrier(state.normalEquationInfo.comm);
    logStatus<<"=== for each receiver separately ============================"<<Log::endl;
    auto estimateSingleReceiver = selectReceivers->select(state.gnss->receivers);
    for(auto recv : state.gnss->receivers)
      if(!recv->useable())
        estimateSingleReceiver.at(recv->idRecv()) = FALSE;
    logInfo<<"  "<<std::count(estimateSingleReceiver.begin(), estimateSingleReceiver.end(), TRUE)<<" receivers selected"<<Log::endl;
    if(Parallel::size(state.normalEquationInfo.comm) > 1)
      logInfo<<"  Only results of a subset of stations are displayed in the following"<<Log::endl;

    // Save old state
    auto normalEquationInfoOld = state.normalEquationInfo;

    // self comm
    state.normalEquationInfo.comm = Parallel::createCommunicator({Parallel::myRank(state.normalEquationInfo.comm)}, state.normalEquationInfo.comm);
    state.normalEquationInfo.isEachReceiverSeparately = TRUE;

    for(UInt idRecv=0; idRecv<state.gnss->receivers.size(); idRecv++)
      if(estimateSingleReceiver.at(idRecv) && state.gnss->receivers.at(idRecv)->isMyRank())
      {
        logStatus<<"=== select single receiver ("<<state.gnss->receivers.at(idRecv)->name()<<") ==========================="<<Log::endl;
        VariableList varList;
        addVariable(variableReceiver, state.gnss->receivers.at(idRecv)->name(), varList);

        std::fill(state.normalEquationInfo.estimateReceiver.begin(), state.normalEquationInfo.estimateReceiver.end(), FALSE);
        state.normalEquationInfo.estimateReceiver.at(idRecv) = TRUE;
        state.changedNormalEquationInfo = TRUE;

        try
        {
          GnssProcessingStepPtr processingSteps;
          configProcessingSteps.read(processingSteps, varList);
          processingSteps->process(state);
        }
        catch(std::exception &e)
        {
          logError<<state.gnss->receivers.at(idRecv)->name()<<": disabled due to exception in single receiver loop:"<<Log::endl;
          logError<<e.what()<<Log::endl;
          state.gnss->receivers.at(idRecv)->disable(e.what());
        }
      } // for(idRecv)

    // restore old state
    std::swap(state.normalEquationInfo, normalEquationInfoOld);
    state.changedNormalEquationInfo = TRUE;

    // synchronize transceivers
    state.gnss->synchronizeTransceivers(state.normalEquationInfo.comm);

    Parallel::barrier(state.normalEquationInfo.comm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
