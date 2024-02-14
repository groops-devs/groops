/***********************************************/
/**
* @file gnssProcessing.cpp
*
* @brief GNSS/LEO satellite orbit determination, station network analysis, PPP.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2010-08-04
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program processes GNSS observations. It calculates the linearized observation equations,
accumulates them into a system of normal equations and solves it.

The primary use cases of this program are:
\begin{itemize}
  \item \reference{GNSS satellite orbit determination and station network analysis}{cookbook.gnssNetwork}
  \item \reference{Kinematic orbit determination of LEO satellites}{cookbook.kinematicOrbit}
  \item \reference{GNSS precise point positioning (PPP)}{cookbook.gnssPpp}
\end{itemize}

The observation epochs are defined by \configClass{timeSeries}{timeSeriesType}
and only observations at these epochs (within a \config{timeMargin}) are considered.

To calculate observation equations from the tracks, the model parameters or unknown parameters need to be
defined beforehand. These unknown parameters can be chosen arbitrarily by the user with an adequate list of defined
\configClass{parametrization}{gnssParametrizationType}.
Some of the \configClass{parametrization}{gnssParametrizationType} also include a priori models.

Lastly it is required to define the process flow of the gnssProcessing. This is accomplished
with a list of \configClass{processingSteps}{gnssProcessingStepType}.
Each step is processed consecutively. Some steps allow the selection of parameters, epochs,
or the normal equation structure, which affects all subsequent steps.
A minimal example consists of following steps:
\begin{itemize}
  \item \configClass{estimate}{gnssProcessingStepType:estimate}: iterative float solution with outlier downeighting
  \item \configClass{resolveAmbiguities}{gnssProcessingStepType:resolveAmbiguities}:
        fix ambiguities to integer and remove them from the normals
  \item \configClass{estimate}{gnssProcessingStepType:estimate}: few iteration for final outlier downweighting
  \item \configClass{writeResults}{gnssProcessingStepType:writeResults}:
        write the output files defined in \configClass{parametrization}{gnssParametrizationType}
\end{itemize}

If the program is run on multiple processes the \configClass{receiver}{gnssReceiverGeneratorType}s
(stations or LEO satellites) are distributed over the processes.

See also \program{GnssSimulateReceiver}.
)";

/***********************************************/

#include "programs/program.h"
#include "config/configRegister.h"
#include "classes/timeSeries/timeSeries.h"
#include "classes/earthRotation/earthRotation.h"
#include "gnss/gnss.h"
#include "gnss/gnssReceiverGenerator/gnssReceiverGenerator.h"
#include "gnss/gnssTransmitterGenerator/gnssTransmitterGenerator.h"
#include "gnss/gnssParametrization/gnssParametrization.h"
#include "gnss/gnssProcessingStep/gnssProcessingStep.h"

/***** CLASS ***********************************/

/** @brief GNSS/LEO satellite orbit determination, station network analysis, PPP.
* @ingroup programsGroup */
class GnssProcessing
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GnssProcessing, PARALLEL, "GNSS/LEO satellite orbit determination, station network analysis, PPP", Gnss)

/***********************************************/

void GnssProcessing::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    TimeSeriesPtr               timeSeries;
    Double                      marginSeconds;
    GnssTransmitterGeneratorPtr transmitterGenerator;
    GnssReceiverGeneratorPtr    receiverGenerator;
    GnssParametrizationPtr      gnssParametrization;
    EarthRotationPtr            earthRotation;
    GnssProcessingStepPtr       processingSteps;

    readConfig(config, "timeSeries",      timeSeries,           Config::MUSTSET,  "",    "defines observation epochs");
    readConfig(config, "timeMargin",      marginSeconds,        Config::DEFAULT,  "0.1", "[seconds] margin to consider two times identical");
    readConfig(config, "transmitter",     transmitterGenerator, Config::MUSTSET,  "",    "constellation of GNSS satellites");
    readConfig(config, "receiver",        receiverGenerator,    Config::MUSTSET,  "",    "ground station network or LEO satellite");
    readConfig(config, "earthRotation",   earthRotation,        Config::MUSTSET,  "",    "apriori earth rotation");
    readConfig(config, "parametrization", gnssParametrization,  Config::MUSTSET,  "",    "models and parameters");
    readConfig(config, "processingStep",  processingSteps,      Config::MUSTSET,  "",    "steps are processed consecutively");
    if(isCreateSchema(config)) return;

    // ============================

    // init the GNSS system
    // --------------------
    logInfo<<"Init GNSS"<<Log::endl;
    std::vector<Time> times = timeSeries->times();
    GnssPtr gnss = std::make_shared<Gnss>();
    gnss->init(times, seconds2time(marginSeconds), transmitterGenerator, receiverGenerator, earthRotation, gnssParametrization, comm);
    receiverGenerator->preprocessing(gnss.get(), comm);
    gnss->synchronizeTransceivers(comm);
    transmitterGenerator = nullptr;
    receiverGenerator    = nullptr;
    gnssParametrization  = nullptr;
    earthRotation        = nullptr;
    logInfo<<"  transmitter: "<<std::count_if(gnss->transmitters.begin(), gnss->transmitters.end(), [](auto t) {return t->useable();})<<Log::endl;
    logInfo<<"  receiver:    "<<std::count_if(gnss->receivers.begin(),    gnss->receivers.end(),    [](auto r) {return r->useable();})<<Log::endl;
    if(!std::any_of(gnss->transmitters.begin(), gnss->transmitters.end(), [](auto trans){return trans->useable();}))
    {
      logWarningOnce<<times.front().dateTimeStr()<<" - "<<times.back().dateTimeStr()<<": no useable transmitters"<<Log::endl;
      return;
    }
    if(!std::any_of(gnss->receivers.begin(), gnss->receivers.end(), [](auto recv){return recv->useable();}))
    {
      logWarningOnce<<times.front().dateTimeStr()<<" - "<<times.back().dateTimeStr()<<": no useable receivers"<<Log::endl;
      return;
    }

    // count observation types
    // -----------------------
    logInfo<<"types and number of observations:"<<Log::endl;
    std::vector<GnssType> types = gnss->types(~(GnssType::PRN + GnssType::FREQ_NO));
    Vector countTypes(types.size());
    for(auto recv : gnss->receivers)
      if(recv->isMyRank())
        for(UInt idEpoch=0; idEpoch<recv->idEpochSize(); idEpoch++)
          for(UInt idTrans=0; idTrans<recv->idTransmitterSize(idEpoch); idTrans++)
          {
            auto obs = recv->observation(idTrans, idEpoch);
            if(obs)
              for(UInt idType=0; idType<obs->size(); idType++)
              {
                const UInt idx = GnssType::index(types, obs->at(idType).type);
                if(idx != NULLINDEX)
                  countTypes(idx)++;
              }
          }
    Parallel::reduceSum(countTypes, 0, comm);

    for(UInt idType=0; idType<types.size(); idType++)
      logInfo<<"  "<<types.at(idType).str()<<":"<<countTypes(idType)%"%10i"s<<Log::endl;
    logInfo<<"        + ========="<<Log::endl;
    logInfo<<"  total:"<<sum(countTypes)%"%11i"s<<Log::endl;

    UInt countTracks = 0;
    for(auto recv : gnss->receivers)
      if(recv->isMyRank())
        countTracks += recv->tracks.size();
    Parallel::reduceSum(countTracks, 0, comm);
    logInfo<<"  number of tracks: "<<countTracks<<Log::endl;

    // Processing steps
    // ----------------
    GnssProcessingStep::State state(gnss, comm);
    processingSteps->process(state);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
