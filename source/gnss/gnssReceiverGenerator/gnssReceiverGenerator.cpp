/***********************************************/
/**
* @file gnssReceiverGenerator.cpp
*
* @brief Provides a list of receivers.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2021-02-25
*
*/
/***********************************************/

#define DOCSTRING_GnssReceiverGenerator

#include "base/import.h"
#include "config/configRegister.h"
#include "config/config.h"
#include "inputOutput/logging.h"
#include "classes/earthRotation/earthRotation.h"
#include "gnss/gnssReceiver.h"
#include "gnss/gnssTransmitter.h"
#include "gnss/gnssReceiverGenerator/gnssReceiverGenerator.h"
#include "gnss/gnssReceiverGenerator/gnssReceiverGeneratorStationNetwork.h"
#include "gnss/gnssReceiverGenerator/gnssReceiverGeneratorLowEarthOrbiter.h"
#include "gnss/gnssReceiverGenerator/gnssReceiverGenerator.h"

/***********************************************/

GROOPS_REGISTER_CLASS(GnssReceiverGenerator, "gnssReceiverGeneratorType",
                      GnssReceiverGeneratorStationNetwork,
                      GnssReceiverGeneratorLowEarthOrbiter)

GROOPS_READCONFIG_UNBOUNDED_CLASS(GnssReceiverGenerator, "gnssReceiverGeneratorType")

/***********************************************/

GnssReceiverGenerator::GnssReceiverGenerator(Config &config, const std::string &name)
{
  try
  {
    std::string type;
    while(readConfigChoice(config, name, type, Config::OPTIONAL, "", ""))
    {
      if(readConfigChoiceElement(config, "stationNetwork", type, ""))
        base.push_back(new GnssReceiverGeneratorStationNetwork(config));
      if(readConfigChoiceElement(config, "lowEarthOrbiter", type, ""))
        base.push_back(new GnssReceiverGeneratorLowEarthOrbiter(config));
      endChoice(config);
      if(isCreateSchema(config))
        return;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

GnssReceiverGenerator::~GnssReceiverGenerator()
{
  for(UInt i=0; i<base.size(); i++)
    delete base.at(i);
}

/***********************************************/

std::vector<GnssReceiverPtr> GnssReceiverGenerator::receivers(std::vector<GnssType> simulationTypes, const std::vector<Time> &times, const Time &timeMargin,
                                                              const std::vector<GnssTransmitterPtr> &transmitters,
                                                              EarthRotationPtr earthRotation, Parallel::CommunicatorPtr comm)
{
  try
  {
    std::vector<GnssReceiverPtr> receivers;
    for(UInt i=0; i<base.size(); i++)
      base.at(i)->init(simulationTypes, times, timeMargin, transmitters, earthRotation, comm, receivers);
    return receivers;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssReceiverGenerator::preprocessing(Gnss *gnss, Parallel::CommunicatorPtr comm)
{
  try
  {
    for(UInt i=0; i<base.size(); i++)
      base.at(i)->preprocessing(gnss, comm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssReceiverGenerator::simulation(NoiseGeneratorPtr noiseClock, NoiseGeneratorPtr noiseObs,
                                       Gnss *gnss, Parallel::CommunicatorPtr comm)
{
  try
  {
    for(UInt i=0; i<base.size(); i++)
      base.at(i)->simulation(noiseClock, noiseObs, gnss, comm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void GnssReceiverGeneratorBase::printPreprocessingInfos(const std::string &header, const std::vector<GnssReceiverPtr> &receivers,
                                                        Bool disabledOnly, Parallel::CommunicatorPtr comm)
{
  try
  {
    // distribute process id of receivers
    // ----------------------------------
    Vector recvProcess(receivers.size());
    for(UInt idRecv=0; idRecv<receivers.size(); idRecv++)
      if(receivers.at(idRecv)->isMyRank() || !receivers.at(idRecv)->disableReason.empty())
        recvProcess(idRecv) = Parallel::myRank(comm)+1;
    Parallel::reduceSum(recvProcess, 0, comm);
    Parallel::broadCast(recvProcess, 0, comm);

    // preprocessing infos
    // -------------------
    for(UInt idRecv=0; idRecv<receivers.size(); idRecv++)
      if(recvProcess(idRecv))
      {
        std::string              reason = receivers.at(idRecv)->disableReason;
        std::vector<std::string> infos  = receivers.at(idRecv)->preprocessingInfos;
        if(recvProcess(idRecv) > 1)
        {
          if(recvProcess(idRecv)-1 == Parallel::myRank(comm))
          {
            Parallel::send(reason, 0, comm);
            Parallel::send(infos,  0, comm);
          }
          if(Parallel::isMaster(comm))
          {
            Parallel::receive(reason, static_cast<UInt>(recvProcess(idRecv)-1), comm);
            Parallel::receive(infos,  static_cast<UInt>(recvProcess(idRecv)-1), comm);
          }
        }

        if(Parallel::isMaster(comm) && !(disabledOnly && reason.empty()))
        {
          logInfo<<receivers.at(idRecv)->name()<<": "<<header<<Log::endl;
          if(!reason.empty())
            logWarning<<receivers.at(idRecv)->name()<<" disabled: "<<reason<<Log::endl;
          for(auto &info : infos)
            logInfo<<"  "<<info<<Log::endl;
        }
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
