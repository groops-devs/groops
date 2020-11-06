/***********************************************/
/**
* @file eclipse.cpp
*
* @brief Shadowing of satellites by moon and Earth.
*
* @author Torsten Mayer-Guerr
* @date 2020-03-08
*
*/
/***********************************************/

#define DOCSTRING_Eclipse

#include "base/import.h"
#include "config/configRegister.h"
#include "classes/eclipse/eclipseConical.h"
#include "classes/eclipse/eclipseSOLAARS.h"
#include "classes/eclipse/eclipse.h"

/***********************************************/

GROOPS_REGISTER_CLASS(Eclipse, "eclipseType",
                      EclipseConical,
                      EclipseSOLAARS)

GROOPS_READCONFIG_CLASS(Eclipse, "eclipseType")

/***********************************************/

EclipsePtr Eclipse::create(Config &config, const std::string &name)
{
  try
  {
    EclipsePtr eclipse;
    std::string choice;
    readConfigChoice(config, name, choice, Config::MUSTSET, "", "Shadowing of satellites by moon and Earth");

    if(readConfigChoiceElement(config, "conical", choice, "Umbra and penumbra shadow model"))
      eclipse = EclipsePtr(new EclipseConical(config));
    if(readConfigChoiceElement(config, "SOLAARS", choice, "Penumbra with Oblateness and Lower Atmospheric Absorption, Refraction, and Scattering"))
      eclipse = EclipsePtr(new EclipseSOLAARS(config));
    endChoice(config);

    return eclipse;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

// For detailed information on this calculation, see Montenbruck and Gill (2000, pp. 80-83)
Double Eclipse::shadowScalingFactor(const Vector3d &posSat, const Vector3d &posSun, const Vector3d &posBody, const Double &radiusBody)
{
  try
  {
    const Double R_Sun       = 6.96342e8;
    const Double distSatSun  = (posSat-posSun).r();
    const Double distSatBody = (posSat-posBody).r();

    // Apparent radius of Sun (a),
    // apparent radius of occulting body (Earth, Moon) (b),
    // apparent separation of the centers of both bodies (c)
    const Double a = std::asin(R_Sun/distSatSun);
    const Double b = std::asin(radiusBody/distSatBody);
    const Double c = std::acos(inner(posBody-posSat, posSun-posSat)/(distSatBody*distSatSun));

    // Satellite is not in occulting body's (Earth, Moon) shadow
    if(a+b <= c)
      return 1.;

    // Satellite is in full occulting body's (Earth, Moon) shadow
    if(std::fabs(a-b) > c)
      return 0.;

    // Satellite is in partial occulting body's (Earth, Moon) shadow
    const Double x = (c*c+a*a-b*b)/(2*c);
    const Double y = std::sqrt(a*a-x*x);
    const Double A = a*a*std::acos(x/a)+b*b*std::acos((c-x)/b)-c*y; // Apparent visible area of the Sun
    return 1-A/(PI*a*a);
  }
  catch (std::exception& e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
