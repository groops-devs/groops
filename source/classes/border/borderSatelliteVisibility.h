/***********************************************/
/**
* @file borderSatelliteVisibility.h
*
* @brief Satellite visibility border.
* @see Border
*
* @author Andreas Kvas
* @date 2026-09-18
*
*/
/***********************************************/

#ifndef __GROOPS_BORDERSATELLITEVISIBILITY__
#define __GROOPS_BORDERSATELLITEVISIBILITY__

// Latex documentation
#ifdef DOCSTRING_Border
static const char *docstringBorderSatelliteVisibility = R"(
\subsection{SatelliteVisibility}

Selects antenna sky directions that are geometrically accessible to a
circular satellite orbit. This border represents the long-term geometrical visibility envelope.
In this context longitude is interpreted as azimuth, measured clockwise from north, and latitude as elevation.

The receiver is represented on a spherical Earth with radius
\config{earthRadius}. The satellite orbit is defined by its height above the Earth's surface
\config{orbitHeight} (added to the Earth's surface) and its \config{inclination}.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "classes/border/border.h"

/***** CLASS ***********************************/

/** @brief Satellite visibility border.
* @ingroup borderGroup
* @see Border */
class BorderSatelliteVisibility : public BorderBase
{
  Angle stationLatitude, minElevation, inclination;
  Double earthRadius, orbitHeight;
  Bool  exclude;

public:
  BorderSatelliteVisibility(Config &config);

  Bool isInnerPoint(Angle lambda, Angle phi) const;
  Bool isExclude() const {return exclude;}
};

/***********************************************/

inline BorderSatelliteVisibility::BorderSatelliteVisibility(Config &config)
{
  readConfig(config, "stationLatitude", stationLatitude, Config::MUSTSET,  "", "geographic latitude of the station in degrees");
  readConfig(config, "earthRadius",     earthRadius,     Config::DEFAULT, STRING_DEFAULT_R, "radius of the spherical Earth in meters");
  readConfig(config, "inclination",     inclination,     Config::MUSTSET,  "", "inclination of the satellite orbit in degrees");
  readConfig(config, "orbitHeight",     orbitHeight,     Config::MUSTSET,  "", "height of the satellite orbit above the Earth's surface in meters");
  readConfig(config, "minElevation",    minElevation,    Config::DEFAULT,  "0", "minimum elevation angle for satellite visibility in degrees");
  readConfig(config, "exclude",         exclude,         Config::DEFAULT,  "0", "dismiss points inside");
  if(isCreateSchema(config)) return;

  Double orbitRadius = earthRadius + orbitHeight;

  if(earthRadius <= 0.0)
    throw(Exception("earthRadius must be positive"));

  if(orbitHeight <= 0.0)
    throw(Exception("orbitHeight must be positive"));

}

/***********************************************/

inline Bool BorderSatelliteVisibility::isInnerPoint(Angle lambda/* interpreted as azimuth*/, Angle phi/* interpreted as elevation*/) const
{
  if(phi < minElevation || phi > Angle(PI/2.0))
    return FALSE;

  const Double R   = earthRadius;
  const Double a   = earthRadius + orbitHeight;
  const Double A   = lambda;
  const Double e   = phi;
  const Double phi0 = stationLatitude;

  const Double sinE = std::sin(e);
  const Double cosE = std::cos(e);

  const Double radicand = a*a - R*R*cosE*cosE;

  if(radicand < 0.0)
    return FALSE;

  const Double distance = -R*sinE + std::sqrt(std::max(0.0, radicand));

  const Double sinBeta = (R*std::sin(phi0) + distance * (sinE*std::sin(phi0) + cosE*std::cos(A)*std::cos(phi0))) / a;

  // abs(sin(i)) also handles retrograde inclinations.
  const Double maximumSinBeta =
      std::abs(std::sin(Double(inclination)));

  return std::abs(sinBeta) <= maximumSinBeta + 1e-14;
}

/***********************************************/

#endif /* __GROOPS_BORDERSATELLITEVISIBILITY__ */
