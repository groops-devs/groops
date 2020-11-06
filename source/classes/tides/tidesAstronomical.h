/***********************************************/
/**
* @file tidesAstronomical.h
*
* @brief Astronomical tides
* Computed from planetary ephemerides.
* @see Tides
*
* @author Torsten Mayer-Guerr
* @date 2015-01-24
*
*/
/***********************************************/

#ifndef __GROOPS_TIDESASTRONOMICAL__
#define __GROOPS_TIDESASTRONOMICAL__

// Latex documentation
#ifdef DOCSTRING_Tides
static const char *docstringTidesAstronomical = R"(
\subsection{AstronomicalTide}\label{tidesType:astronomicalTide}
This class computes the tide generating potential (TGP) of sun, moon
and planets (Mercury, Venus, Mars, Jupiter, Saturn).
It takes into account the flattening of the Earth (At the moment only at the acceleration level).

The computed result is multiplied with \config{factor}.
)";
#endif


/***********************************************/

#include "classes/tides/tides.h"

/***** CLASS ***********************************/

/** @brief Astronomical tides.
* @ingroup tidesGroup
* Computed from planetary ephemerids.
* @see Tides */
class TidesAstronomical : public TidesBase
{
  Double             factor;
  Bool               useEarth, useMoon, useSun, usePlanets;
  SphericalHarmonics j2earth;

  Double   directTidePotential     (Double GM, const Vector3d &point1, const Vector3d &point2) const;
  Double   directTideRadialGradient(Double GM, const Vector3d &point1, const Vector3d &point2) const;
  Vector3d directTideAcceleration  (Double GM, const Vector3d &point1, const Vector3d &point2) const;
  Tensor3d directTideGradient      (Double GM, const Vector3d &point1, const Vector3d &point2) const;

public:
  TidesAstronomical(Config &config);

  Double   potential      (const Time &time, const Vector3d &point, const Rotary3d &rotEarth, EarthRotationPtr rotation, EphemeridesPtr ephemerides) const override;
  Double   radialGradient (const Time &time, const Vector3d &point, const Rotary3d &rotEarth, EarthRotationPtr rotation, EphemeridesPtr ephemerides) const override;
  Vector3d gravity        (const Time &time, const Vector3d &point, const Rotary3d &rotEarth, EarthRotationPtr rotation, EphemeridesPtr ephemerides) const override;
  Tensor3d gravityGradient(const Time &time, const Vector3d &point, const Rotary3d &rotEarth, EarthRotationPtr rotation, EphemeridesPtr ephemerides) const override;
  Vector3d deformation    (const Time &time, const Vector3d &point, const Rotary3d &rotEarth, EarthRotationPtr rotation, EphemeridesPtr ephemerides,
                           Double gravity, const Vector &hn, const Vector &ln) const override;
  void     deformation    (const std::vector<Time> &time, const std::vector<Vector3d> &point, const std::vector<Rotary3d> &rotEarth,
                           EarthRotationPtr rotation, EphemeridesPtr ephemerides, const std::vector<Double> &gravity, const Vector &hn, const Vector &ln,
                           std::vector<std::vector<Vector3d>> &disp) const override;
  SphericalHarmonics sphericalHarmonics(const Time &time, const Rotary3d &rotEarth, EarthRotationPtr rotation, EphemeridesPtr ephemerides,
                                        UInt maxDegree, UInt minDegree, Double GM, Double R) const override;
};

/***********************************************/

#endif
