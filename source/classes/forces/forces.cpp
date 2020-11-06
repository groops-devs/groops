/***********************************************/
/**
* @file forces.cpp
*
* @brief Wraps forces and force-generating potentials.
*
* @author Matthias Ellmer
* @date 2017-01-18
*
*/
/***********************************************/

#define DOCSTRING_Forces

#include "base/import.h"
#include "config/configRegister.h"
#include "classes/forces/forces.h"

/***********************************************/

GROOPS_REGISTER_CLASS_WITHOUT_SUBS(Forces, "forcesType")
GROOPS_READCONFIG_CLASS(Forces, "forcesType")

/***********************************************/

Forces::Forces(Config &config, const std::string &name)
{
  try
  {
    readConfigSequence(config, name, Config::MUSTSET, "", "");
    readConfig(config, "gravityfield",      gravityfield,      Config::OPTIONAL, "", "");
    readConfig(config, "tides",             tides,             Config::OPTIONAL, "", "");
    readConfig(config, "miscAccelerations", miscAccelerations, Config::OPTIONAL, "", "");
    endSequence(config);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector3d Forces::acceleration(SatelliteModelPtr satellite, const Time &time, const Vector3d &position, const Vector3d &velocity,
                              const Rotary3d &rotSat, const Rotary3d &rotEarth, EarthRotationPtr rotation, EphemeridesPtr ephemerides) const
{
  try
  {
    const Vector3d posEarth = rotEarth.rotate(position);
    Vector3d g;
    if(gravityfield)      g += gravityfield->gravity(time, posEarth);
    if(tides)             g += tides->acceleration(time, posEarth, rotEarth, rotation, ephemerides);
    if(miscAccelerations) g += miscAccelerations->acceleration(satellite, time, position, velocity, rotSat, rotEarth, ephemerides);
    return g;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
