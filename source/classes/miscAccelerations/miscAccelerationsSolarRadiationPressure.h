/***********************************************/
/**
* @file miscAccelerationsSolarRadiationPressure.h
*
* @brief DEPRECATED. Use radiationPressure instead.
* @see MiscAccelerations
*
* @author Torsten Mayer-Guerr
* @date 2014-21-10
*
*/
/***********************************************/

#ifndef __GROOPS_MISCACCELERATIONSSOLARRADIATIONPRESSURE__
#define __GROOPS_MISCACCELERATIONSSOLARRADIATIONPRESSURE__

// Latex documentation
#ifdef DOCSTRING_MiscAccelerations
static const char *docstringMiscAccelerationsSolarRadiationPressure = R"(
\subsection{SolarRadiationPressure}\label{miscAccelerationsType:solarRadiationPressure}
DEPRECATED. Use radiationPressure instead.
)";
#endif

/***********************************************/

#include "classes/eclipse/eclipse.h"
#include "classes/miscAccelerations/miscAccelerationsRadiationPressure.h"
#include "classes/miscAccelerations/miscAccelerations.h"

/***** CLASS ***********************************/

/** @brief DEPRECATED. Use radiationPressure instead.
 * @ingroup miscAccelerationsGroup
 * @see MiscAccelerations */
class MiscAccelerationsSolarRadiationPressure : public MiscAccelerationsBase
{
  EclipsePtr eclipse;
  Double     solarflux;
  Double     factor;

public:
  MiscAccelerationsSolarRadiationPressure(Config &config);

  Vector3d acceleration(SatelliteModelPtr satellite, const Time &time, const Vector3d &position, const Vector3d &velocity,
                        const Rotary3d &rotSat, const Rotary3d &rotEarth, EphemeridesPtr ephemerides) override;
};

/***********************************************/

inline MiscAccelerationsSolarRadiationPressure::MiscAccelerationsSolarRadiationPressure(Config &config)
{
  try
  {
    readConfig(config, "solarflux", solarflux, Config::DEFAULT, "1367", "solar flux constant in 1 AU [W/m**2]");
    readConfig(config, "eclipse",   eclipse,   Config::MUSTSET, "",     "");
    readConfig(config, "factor",    factor,    Config::DEFAULT, "1.0",  "the result is multiplied by this factor, set -1 to subtract the field");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Vector3d MiscAccelerationsSolarRadiationPressure::acceleration(SatelliteModelPtr satellite, const Time &time,
                                                                      const Vector3d &position, const Vector3d &/*velocity*/,
                                                                      const Rotary3d &rotSat, const Rotary3d &rotEarth, EphemeridesPtr ephemerides)
{
  try
  {
    if(!satellite)
      throw(Exception("No satellite model given"));
    if(!ephemerides)
      throw(Exception("No ephemerides given"));

    // from sun to satellite in SRF
    Vector3d direction = rotSat.inverseRotate(position - ephemerides->position(time, Ephemerides::SUN));
    const Double distanceSun = direction.normalize();
    const Double AU = 149597870700.0;
    const Double ny = eclipse ? eclipse->factor(time, position, ephemerides) : 1; // shadowing from Earth and Moon

    Vector absorbedPressure(satellite->surfaces.size()); // power for each satellite surface
    Vector3d F = MiscAccelerationsRadiationPressure::force(satellite, direction, ny*solarflux/LIGHT_VELOCITY*std::pow(AU/distanceSun,2), 0., absorbedPressure);

    // compute thermal radition
    for(UInt i=0; i<satellite->surfaces.size(); i++)
    {
      auto &surface = satellite->surfaces.at(i);
      if(surface.specificHeatCapacity < 0)
        F -= (2./3.) * surface.area * absorbedPressure(i) * surface.normal; // spontaneous reemission of absorbed radiation
    }

    return (factor/satellite->mass) * rotEarth.rotate(rotSat.rotate(F));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
