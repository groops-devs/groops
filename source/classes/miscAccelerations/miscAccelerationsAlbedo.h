/***********************************************/
/**
* @file miscAccelerationsAlbedo.h
*
* @brief DEPRECATED. Use radiationPressure instead.
* @see MiscAccelerations
*
* @author Torsten Mayer-Guerr
* @date 2014-21-10
*
*/
/***********************************************/

#ifndef __GROOPS_MISCACCELERATIONALBEDO__
#define __GROOPS_MISCACCELERATIONALBEDO__

// Latex documentation
#ifdef DOCSTRING_MiscAccelerations
static const char *docstringMiscAccelerationsAlbedo = R"(
\subsection{Albedo}\label{miscAccelerationsType:albedo}
DEPRECATED. Use radiationPressure instead.
)";
#endif

/***********************************************/

#include "base/planets.h"
#include "classes/miscAccelerations/miscAccelerations.h"

/***** CLASS ***********************************/

/** @brief DEPRECATED. Use radiationPressure instead.
* @ingroup miscAccelerationsGroup
* @see MiscAccelerations */
class MiscAccelerationsAlbedo : public MiscAccelerationsBase
{
  std::vector<Vector3d>            points;
  std::vector<Double>              areas;
  std::vector<std::vector<Double>> reflectivity;
  std::vector<std::vector<Double>> emissivity;
  Double                           solarflux;
  Double                           factor;

public:
  MiscAccelerationsAlbedo(Config &config);

  Vector3d acceleration(SatelliteModelPtr satellite, const Time &time, const Vector3d &position, const Vector3d &velocity,
                        const Rotary3d &rotSat, const Rotary3d &rotEarth, EphemeridesPtr ephemerides) override;
};

/***********************************************/

#endif
