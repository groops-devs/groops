/***********************************************/
/**
* @file magnetosphere.cpp
*
* @brief Magentic field of the Earth.
*
* @author Torsten Mayer-Guerr
* @date 2020-06-08
*
*/
/***********************************************/

#define DOCSTRING_Magnetosphere

#include "base/import.h"
#include "base/planets.h"
#include "config/configRegister.h"
#include "classes/magnetosphere/magnetosphereIgrf.h"
#include "classes/magnetosphere/magnetosphere.h"

/***********************************************/

GROOPS_REGISTER_CLASS(Magnetosphere, "magnetosphereType",
                      MagnetosphereIgrf)
GROOPS_READCONFIG_CLASS(Magnetosphere, "magnetosphereType")

/***********************************************/

MagnetospherePtr Magnetosphere::create(Config &config, const std::string &name)
{
  try
  {
    MagnetospherePtr ptr;
    std::string choice;

    readConfigChoice(config, name, choice, Config::MUSTSET, "", "magentic field of the Earth");
    if(readConfigChoiceElement(config, "igrf", choice, "International Geomagnetic Reference Field"))
      ptr = std::make_shared<MagnetosphereIgrf>(config);
    endChoice(config);

    return ptr;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Rotary3d Magnetosphere::rotaryCelestial2SolarGeomagneticFrame(const Time &time) const
{
  try
  {
    const Vector3d z = Planets::celestial2TerrestrialFrame(time).inverseRotate(geomagneticNorthPole(time));
    const Vector3d y = normalize(crossProduct(z, Planets::positionSun(time)));
    const Vector3d x = crossProduct(y, z);
    return inverse(Rotary3d(x,y));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
