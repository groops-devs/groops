/***********************************************/
/**
* @file gnssParametrizationReceiver.cpp
*
* @brief GNSS receivers.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2020-06-02
*
*/
/***********************************************/

#define DOCSTRING_GnssParametrizationReceiver

#include "config/configRegister.h"
#include "gnss/gnss.h"
#include "gnss/gnssParametrizationReceiverStationNetwork.h"
#include "gnss/gnssParametrizationReceiverLeo.h"
#include "gnss/gnssParametrizationReceiver.h"

/***********************************************/

GROOPS_REGISTER_CLASS(GnssParametrizationReceiver, "gnssParametrizationReceiverType",
                      GnssParametrizationReceiverStationNetwork,
                      GnssParametrizationReceiverLeo)

GROOPS_READCONFIG_CLASS(GnssParametrizationReceiver, "gnssParametrizationReceiverType")

/***********************************************/

GnssParametrizationReceiverPtr GnssParametrizationReceiver::create(Config &config, const std::string &name)
{
  try
  {
    GnssParametrizationReceiverPtr ptr;
    std::string      type;

    readConfigChoice(config, name, type, Config::MUSTSET, "", "");
    if(readConfigChoiceElement(config, "stationNetwork", type, ""))
      ptr = std::make_shared<GnssParametrizationReceiverStationNetwork>(config);
    if(readConfigChoiceElement(config, "lowEarthOrbiter", type, ""))
      ptr = std::make_shared<GnssParametrizationReceiverLeo>(config);
    endChoice(config);
    return ptr;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
