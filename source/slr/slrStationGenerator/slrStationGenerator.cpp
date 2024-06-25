/***********************************************/
/**
* @file slrStationGenerator.cpp
*
* @brief Provides a list of stations.
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#define DOCSTRING_SlrStationGenerator

#include "base/import.h"
#include "config/configRegister.h"
#include "config/config.h"
#include "inputOutput/logging.h"
#include "classes/earthRotation/earthRotation.h"
#include "slr/slrStation.h"
#include "slr/slrSatellite.h"
#include "slr/slrStationGenerator/slrStationGenerator.h"
#include "slr/slrStationGenerator/slrStationGeneratorStations.h"
#include "slr/slrStationGenerator/slrStationGenerator.h"

/***********************************************/

GROOPS_REGISTER_CLASS(SlrStationGenerator, "slrStationGeneratorType",
                      SlrStationGeneratorStations)

GROOPS_READCONFIG_UNBOUNDED_CLASS(SlrStationGenerator, "slrStationGeneratorType")

/***********************************************/

SlrStationGenerator::SlrStationGenerator(Config &config, const std::string &name)
{
  try
  {
    std::string type;
    while(readConfigChoice(config, name, type, Config::OPTIONAL, "", ""))
    {
      if(readConfigChoiceElement(config, "stations", type, ""))
        base.push_back(new SlrStationGeneratorStations(config));
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

SlrStationGenerator::~SlrStationGenerator()
{
  for(UInt i=0; i<base.size(); i++)
    delete base.at(i);
}

/***********************************************/

std::vector<SlrStationPtr> SlrStationGenerator::stations(const std::vector<Time> &times, const std::vector<SlrSatellitePtr> &satellites,
                                                         EarthRotationPtr earthRotation)
{
  try
  {
    std::vector<SlrStationPtr> stations;
    for(UInt i=0; i<base.size(); i++)
      base.at(i)->init(times, satellites, earthRotation, stations);
    return stations;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

// void SlrStationGenerator::preprocessing(Slr *slr)
// {
//   try
//   {
//     for(UInt i=0; i<base.size(); i++)
//       base.at(i)->preprocessing(slr);
//   }
//   catch(std::exception &e)
//   {
//     GROOPS_RETHROW(e)
//   }
// }

/***********************************************/

// void SlrStationGenerator::simulation(const std::vector<SlrType> &types, NoiseGeneratorPtr noiseClock, NoiseGeneratorPtr noiseObs,
//                                        Slr *slr)
// {
//   try
//   {
//     for(UInt i=0; i<base.size(); i++)
//       base.at(i)->simulation(types, noiseClock, noiseObs, slr);
//   }
//   catch(std::exception &e)
//   {
//     GROOPS_RETHROW(e)
//   }
// }

/***********************************************/
