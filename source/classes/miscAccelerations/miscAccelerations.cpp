/***********************************************/
/**
* @file miscAccelerations.cpp
*
* @brief Non conservative forces acting on satellites.
*
* @author Torsten Mayer-Guerr
* @date 2013-02-09
*
*/
/***********************************************/

#define DOCSTRING_MiscAccelerations

#include "base/import.h"
#include "config/configRegister.h"
#include "classes/miscAccelerations/miscAccelerationsRelativisticEffect.h"
#include "classes/miscAccelerations/miscAccelerationsSolarRadiationPressure.h"
#include "classes/miscAccelerations/miscAccelerationsAtmosphericDrag.h"
#include "classes/miscAccelerations/miscAccelerationsAlbedo.h"
#include "classes/miscAccelerations/miscAccelerationsAntennaThrust.h"
#include "classes/miscAccelerations/miscAccelerationsFromParametrization.h"
#include "classes/miscAccelerations/miscAccelerations.h"

/***********************************************/

GROOPS_REGISTER_CLASS(MiscAccelerations, "miscAccelerationsType",
                      MiscAccelerationsRelativisticEffect,
                      MiscAccelerationsSolarRadiationPressure,
                      MiscAccelerationsAlbedo,
                      MiscAccelerationsAtmosphericDrag,
                      MiscAccelerationsAntennaThrust,
                      MiscAccelerationsFromParametrization)

GROOPS_READCONFIG_UNBOUNDED_CLASS(MiscAccelerations, "miscAccelerationsType")

/***********************************************/

MiscAccelerations::MiscAccelerations(Config &config, const std::string &name)
{
  try
  {
    std::string type;
    while(readConfigChoice(config, name, type, Config::OPTIONAL, "", "non conservative forces acting on satellites"))
    {
      if(readConfigChoiceElement(config, "relativisticEffect", type, "Relativistic effect (IERS2010)"))
       acc.push_back(new MiscAccelerationsRelativisticEffect(config));
      if(readConfigChoiceElement(config, "solarRadiationPressure", type, "Solar radiation pressure model"))
       acc.push_back(new MiscAccelerationsSolarRadiationPressure(config));
      if(readConfigChoiceElement(config, "albedo", type, "Albedo radiation"))
       acc.push_back(new MiscAccelerationsAlbedo(config));
      if(readConfigChoiceElement(config, "atmosphericDrag", type, "atmospheric drag"))
       acc.push_back(new MiscAccelerationsAtmosphericDrag(config));
      if(readConfigChoiceElement(config, "antennaThrust", type, "antenna thrust"))
       acc.push_back(new MiscAccelerationsAntennaThrust(config));
      if(readConfigChoiceElement(config, "fromParametrization", type, "from a solution vector with given parametrization"))
       acc.push_back(new MiscAccelerationsFromParametrization(config));
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

MiscAccelerations::~MiscAccelerations()
{
  for(UInt i=0; i<acc.size(); i++)
    delete acc.at(i);
}

/***********************************************/

Vector3d MiscAccelerations::acceleration(SatelliteModelPtr satellite, const Time &time,
                                         const Vector3d &position, const Vector3d &velocity,
                                         const Rotary3d &rotSat, const Rotary3d &rotEarth, EphemeridesPtr ephemerides)
{
  try
  {
    if(satellite)
    {
      Vector3d positionSun;
      if(ephemerides)
        positionSun = ephemerides->position(time, Ephemerides::SUN);
      satellite->changeState(time, position, velocity, positionSun, rotSat, rotEarth);
    }

    Vector3d sum;
    for(UInt i=0; i<acc.size(); i++)
      sum += acc.at(i)->acceleration(satellite, time, position, velocity, rotSat, rotEarth, ephemerides);
    return sum;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
