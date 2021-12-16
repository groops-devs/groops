/***********************************************/
/**
* @file gnssSimulateReceiver.cpp
*
* @brief GNSS receiver simulation.
*
* @author Sebastian Strasser
* @author Torsten Mayer-Gürr
* @date 2014-09-25
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program simulates observations from receivers to GNSS satellites.
These simulated observations can then be used in \program{GnssProcessing}, for example to conduct closed-loop simulations.

One or more GNSS constellations must be defined via \configClass{transmitter}{gnssTransmitterGeneratorType}.
Receivers such as ground station networks or Low Earth Orbit (LEO) satellites can be defined via \configClass{receiver}{gnssReceiverGeneratorType}.

If multiple receivers defined an \configFile{outputfileGnssReceiver}{instrument} and \configFile{outputfileClock}{instrument}
are written for each single receiver with the \reference{variable}{general.parser} \verb|{station}| being replaced by the receiver name.

A list of simulated observation types can be defined via \configClass{observationType}{gnssType}. Noise can be added to both observations and clock errors
via \configClass{noiseObervation}{noiseGeneratorType} and \configClass{noiseClockReceiver}{noiseGeneratorType}, respectively. Observation noise is
interpreted as a factor that is multiplied to the accuracy derived from the accuracy pattern of the respective observation type
(see \configFile{inputfileAccuracyDefinition}{gnssAntennaDefinition} in \configClass{receiver}{gnssReceiverGeneratorType}).

The \configClass{parametrization}{gnssParametrizationType} are used to simulate a priori models (e.g. troposphere, signal biases).
Parameter settings and outputfiles are irgnored.

If the program is run on multiple processes the \configClass{receiver}{gnssReceiverGeneratorType}s
(stations or LEO satellites) are distributed over the processes.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/noiseGenerator/noiseGenerator.h"
#include "classes/timeSeries/timeSeries.h"
#include "gnss/gnss.h"
#include "gnss/gnssObservation.h"
#include "gnss/gnssTransmitterGenerator/gnssTransmitterGenerator.h"
#include "gnss/gnssReceiverGenerator/gnssReceiverGenerator.h"

/***** CLASS ***********************************/

/** @brief GNSS receiver simulation.
* @ingroup programsGroup */
class GnssSimulateReceiver
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GnssSimulateReceiver, PARALLEL, "GNSS receiver simulation", Gnss, Simulation)

/***********************************************/

void GnssSimulateReceiver::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName                    fileNameReceiver, fileNameClock;
    TimeSeriesPtr               timeSeries;
    Double                      marginSeconds;
    GnssTransmitterGeneratorPtr transmitterGenerator;
    GnssReceiverGeneratorPtr    receiverGenerator;
    GnssParametrizationPtr      gnssParametrization;
    EarthRotationPtr            earthRotation;
    NoiseGeneratorPtr           noiseClock, noiseObs;
    std::vector<GnssType>       obsTypes;

    readConfig(config, "outputfileGnssReceiver",  fileNameReceiver,     Config::MUSTSET,  "gnssReceiver_{loopTime:%D}.{station}.dat.gz", "variable {station} available, simulated observations");
    readConfig(config, "outputfileClock",         fileNameClock,        Config::OPTIONAL, "clock_{loopTime:%D}.{station}.dat", "variable {station} available, simulated receiver clock errors");
    readConfig(config, "timeSeries",              timeSeries,           Config::MUSTSET,  "",    "defines observation epochs");
    readConfig(config, "timeMargin",              marginSeconds,        Config::DEFAULT,  "0.1", "[seconds] margin to consider two times identical");
    readConfig(config, "transmitter",             transmitterGenerator, Config::MUSTSET,  "",    "constellation of GNSS satellites");
    readConfig(config, "receiver",                receiverGenerator,    Config::MUSTSET,  "",    "ground station network or LEO satellite");
    readConfig(config, "earthRotation",           earthRotation,        Config::MUSTSET,  "",    "apriori earth rotation");
    readConfig(config, "parametrization",         gnssParametrization,  Config::DEFAULT,  "1",   "models and parameters");
    readConfig(config, "observationType",         obsTypes,             Config::MUSTSET,  "",    "simulated observation types");
    readConfig(config, "noiseObservation",        noiseObs,             Config::DEFAULT,  "",    "[-] noise is multiplied with type accuracy pattern of receiver");
    readConfig(config, "noiseClockReceiver",      noiseClock,           Config::DEFAULT,  "",    "[m] noise added to the simulated receiver clock");
    if(isCreateSchema(config)) return;

    // ============================

    // init the GNSS system
    // --------------------
    logInfo<<"Init GNSS"<<Log::endl;
    std::sort(obsTypes.begin(), obsTypes.end());
    std::vector<Time> times = timeSeries->times();
    Gnss gnss;
    gnss.init(times, seconds2time(marginSeconds), transmitterGenerator, receiverGenerator, earthRotation, gnssParametrization, comm);
    receiverGenerator->simulation(obsTypes, noiseClock, noiseObs, &gnss, comm);
    gnss.synchronizeTransceivers(comm);
    logInfo<<"  transmitter: "<<std::count_if(gnss.transmitters.begin(), gnss.transmitters.end(), [](auto t) {return t->useable();})<<Log::endl;
    logInfo<<"  receiver:    "<<std::count_if(gnss.receivers.begin(),    gnss.receivers.end(),    [](auto r) {return r->useable();})<<Log::endl;
    if(!std::any_of(gnss.transmitters.begin(), gnss.transmitters.end(), [](auto trans){return trans->useable();}))
    {
      logWarningOnce<<times.front().dateTimeStr()<<" - "<<times.back().dateTimeStr()<<": no useable transmitters"<<Log::endl;
      return;
    }
    if(!std::any_of(gnss.receivers.begin(), gnss.receivers.end(), [](auto recv){return recv->useable();}))
    {
      logWarningOnce<<times.front().dateTimeStr()<<" - "<<times.back().dateTimeStr()<<": no useable receivers"<<Log::endl;
      return;
    }

    // count observation types
    // -----------------------
    logInfo<<"types and number of observations:"<<Log::endl;
    std::vector<GnssType> types = gnss.types(~(GnssType::PRN + GnssType::FREQ_NO));
    Vector countTypes(types.size());
    for(auto recv : gnss.receivers)
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
    for(auto recv : gnss.receivers)
      if(recv->isMyRank())
        countTracks += recv->tracks.size();
    Parallel::reduceSum(countTracks, 0, comm);
    logInfo<<"  number of tracks: "<<countTracks<<Log::endl;

    // ============================

    // Write observations
    // ------------------
    if(!fileNameReceiver.empty())
    {
      VariableList fileNameVariableList;
      addVariable("station", "****", fileNameVariableList);
      logStatus<<"write receiver observations to files <"<<fileNameReceiver(fileNameVariableList)<<">"<<Log::endl;
      for(auto recv : gnss.receivers)
        if(recv->isMyRank())
        {
          GnssReceiverArc arc;
          for(UInt idEpoch=0; idEpoch<gnss.times.size(); idEpoch++)
            if(recv->useable(idEpoch))
            {
              GnssReceiverEpoch epoch;
              epoch.time       = gnss.times.at(idEpoch);
              epoch.clockError = recv->clockError(idEpoch);

              // get types
              for(UInt idTrans=0; idTrans<recv->idTransmitterSize(idEpoch); idTrans++)
                if(recv->observation(idTrans, idEpoch) && gnss.transmitters.at(idTrans)->useable(idEpoch))
                  for(UInt idType=0; idType<recv->observation(idTrans, idEpoch)->size(); idType++)
                    if(!recv->observation(idTrans, idEpoch)->at(idType).type.isInList(epoch.obsType))
                      epoch.obsType.push_back(recv->observation(idTrans, idEpoch)->at(idType).type & ~(GnssType::PRN+GnssType::FREQ_NO));
              std::sort(epoch.obsType.begin(), epoch.obsType.end());
              if(!epoch.obsType.size())
                continue;

              for(UInt idTrans=0; idTrans<recv->idTransmitterSize(idEpoch); idTrans++)
                if(recv->observation(idTrans, idEpoch) && gnss.transmitters.at(idTrans)->useable(idEpoch))
                {
                  const GnssObservation &obs = *recv->observation(idTrans, idEpoch);
                  const GnssType prn = obs.at(0).type & (GnssType::SYSTEM + GnssType::PRN + GnssType::FREQ_NO);
                  UInt idType = std::distance(epoch.obsType.begin(), std::find(epoch.obsType.begin(), epoch.obsType.end(), prn));

                  epoch.satellite.push_back(prn);
                  for(; (idType<epoch.obsType.size()) && (epoch.obsType.at(idType) == prn); idType++)
                  {
                    epoch.observation.push_back(NAN_EXPR);
                    for(UInt i=0; i<obs.size(); i++)
                      if(obs.at(i).type == epoch.obsType.at(idType))
                      {
                        epoch.observation.back() = obs.at(i).observation;
                        break;
                      }
                  }
                } // for(idTrans)

              if(epoch.satellite.size())
                arc.push_back(epoch);
            } // for(idEpoch)

          fileNameVariableList["station"]->setValue(recv->name());
          InstrumentFile::write(fileNameReceiver(fileNameVariableList), arc);
        } // for(recv)
    } // if(fileNameReceiver)

    // ============================

    // write clock errors
    // ------------------
    if(!fileNameClock.empty())
    {
      VariableList fileNameVariableList;
      addVariable("station", "****", fileNameVariableList);
      logStatus<<"write receiver clocks to files <"<<fileNameClock(fileNameVariableList)<<">"<<Log::endl;
      for(auto recv : gnss.receivers)
        if(recv->isMyRank())
        {
          MiscValueArc arc;
          for(UInt idEpoch=0; idEpoch<gnss.times.size(); idEpoch++)
            if(recv->useable(idEpoch))
            {
              MiscValueEpoch epoch;
              epoch.time  = gnss.times.at(idEpoch);
              epoch.value = recv->clockError(idEpoch);
              arc.push_back(epoch);
            }
          fileNameVariableList["station"]->setValue(recv->name());
          InstrumentFile::write(fileNameClock(fileNameVariableList), arc);
        } // for(recv)
    } // if(fileNameClock)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
