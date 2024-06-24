/***********************************************/
/**
* @file slrSatelliteGenerator.cpp
*
* @brief Provides a list of satellites.
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#define DOCSTRING_SlrSatelliteGenerator

#include "base/import.h"
#include "config/configRegister.h"
#include "config/config.h"
#include "slr/slrSatellite.h"
#include "slr/slrSatelliteGenerator/slrSatelliteGenerator.h"
#include "slr/slrSatelliteGenerator/slrSatelliteGeneratorSatellites.h"

/***********************************************/

GROOPS_REGISTER_CLASS(SlrSatelliteGenerator, "slrSatelliteGeneratorType",
                      SlrSatelliteGeneratorSatellites)

GROOPS_READCONFIG_UNBOUNDED_CLASS(SlrSatelliteGenerator, "slrSatelliteGeneratorType")

/***********************************************/

SlrSatelliteGenerator::SlrSatelliteGenerator(Config &config, const std::string &name)
{
  try
  {
    std::string type;
    while(readConfigChoice(config, name, type, Config::OPTIONAL, "", ""))
    {
      if(readConfigChoiceElement(config, "satellites", type, ""))
        base.push_back(new SlrSatelliteGeneratorSatellites(config));
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

SlrSatelliteGenerator::~SlrSatelliteGenerator()
{
  for(UInt i=0; i<base.size(); i++)
    delete base.at(i);
}

/***********************************************/

std::vector<SlrSatellitePtr> SlrSatelliteGenerator::satellites(const std::vector<Time> &times)
{
  try
  {
    std::vector<SlrSatellitePtr> satellites;
    for(UInt i=0; i<base.size(); i++)
      base.at(i)->init(times, satellites);
    return satellites;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
