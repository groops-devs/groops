/***********************************************/
/**
* @file gnssTransmitterGenerator.cpp
*
* @brief Provides a list of transmitters.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2021-02-25
*
*/
/***********************************************/

#define DOCSTRING_GnssTransmitterGenerator

#include "base/import.h"
#include "config/configRegister.h"
#include "config/config.h"
#include "gnss/gnssTransmitter.h"
#include "gnss/gnssTransmitterGenerator/gnssTransmitterGenerator.h"
#include "gnss/gnssTransmitterGenerator/gnssTransmitterGeneratorGnss.h"

/***********************************************/

GROOPS_REGISTER_CLASS(GnssTransmitterGenerator, "gnssTransmitterGeneratorType",
                      GnssTransmitterGeneratorGnss)

GROOPS_READCONFIG_UNBOUNDED_CLASS(GnssTransmitterGenerator, "gnssTransmitterGeneratorType")

/***********************************************/

GnssTransmitterGenerator::GnssTransmitterGenerator(Config &config, const std::string &name)
{
  try
  {
    std::string type;
    while(readConfigChoice(config, name, type, Config::OPTIONAL, "", ""))
    {
      if(readConfigChoiceElement(config, "GNSS", type, ""))
        base.push_back(new GnssTransmitterGeneratorGnss(config));
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

GnssTransmitterGenerator::~GnssTransmitterGenerator()
{
  for(UInt i=0; i<base.size(); i++)
    delete base.at(i);
}

/***********************************************/

std::vector<GnssTransmitterPtr> GnssTransmitterGenerator::transmitters(const std::vector<Time> &times)
{
  try
  {
    std::vector<GnssTransmitterPtr> transmitters;
    for(UInt i=0; i<base.size(); i++)
      base.at(i)->init(times, transmitters);
    return transmitters;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
