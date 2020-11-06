/***********************************************/
/**
* @file miscAccelerationsAlbedo.h
*
* @brief Albedo radiation.
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
Acceleration caused by Earth's albedo.

The acceleration on the satellite caused by reflection and reradiation is computed using the
algorithm in:

Knocke, P. C., Ries, J. C., and Tapley, B. D. (1988). Earth radiation pressure effects on satellites.
Proceedings of the AIAA/AAS Astrodynamics Conference, 88-4292-CP, 577-87. DOI: 10.2514/6.
1988-4292.
)";
#endif

/***********************************************/

#include "base/planets.h"
#include "classes/miscAccelerations/miscAccelerations.h"

/***** CLASS ***********************************/

/** @brief Albedo radiation.
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
