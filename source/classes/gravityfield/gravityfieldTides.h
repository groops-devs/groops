/***********************************************/
/**
* @file gravityfieldTides.h
*
* @brief Tidal forces computed with class Tide.
* @see Gravityfield
* @see Tides
*
* @author Torsten Mayer-Guerr
* @date 2004-11-01
*
*/
/***********************************************/

#ifndef __GROOPS_GRAVITYFIELDTIDES__
#define __GROOPS_GRAVITYFIELDTIDES__

// Latex documentation
#ifdef DOCSTRING_Gravityfield
static const char *docstringGravityfieldTides = R"(
\subsection{Tides}\label{gravityfieldType:tides}
Treat \configClass{tides}{tidesType} as gravitational forces.
The tides need a realization of \configClass{earthRotation}{earthRotationType}
to transform between the CRF and TRF and to compute rotational deformation
from polar motion.
It also needs \configClass{ephemerides}{ephemeridesType} from Sun, moon, and planets.
)";
#endif

/***********************************************/

#include "classes/tides/tides.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/ephemerides/ephemerides.h"
#include "classes/gravityfield/gravityfield.h"

/***** CLASS ***********************************/

/** @brief Tidal forces computed with class Tide.
* @ingroup gravityfieldGroup
* @see Gravityfield
* @see Tides */
class GravityfieldTides : public GravityfieldBase
{
  EphemeridesPtr   ephemerides;
  EarthRotationPtr earthRotation;
  TidesPtr         tides;

public:
  GravityfieldTides(Config &config);

  Double   potential      (const Time &time, const Vector3d &point) const override;
  Double   radialGradient (const Time &time, const Vector3d &point) const override;
  Vector3d gravity        (const Time &time, const Vector3d &point) const override;
  Tensor3d gravityGradient(const Time &time, const Vector3d &point) const override;
  Vector3d deformation    (const Time &time, const Vector3d &point, Double gravity, const Vector &hn, const Vector &ln) const override;
  void     deformation    (const std::vector<Time> &time, const std::vector<Vector3d> &point, const std::vector<Double> &gravity,
                           const Vector &hn, const Vector &ln, std::vector<std::vector<Vector3d>> &disp) const override;

  SphericalHarmonics sphericalHarmonics(const Time &time, UInt maxDegree, UInt minDegree, Double GM, Double R) const override;

  void variance(const Time &time, const std::vector<Vector3d> &point, const Kernel &kernel, Matrix &D) const override;
};

/***********************************************/

#endif /* __GROOPS_GRAVITYFIELD__ */
