/***********************************************/
/**
 * @file tidesGroup.h
 *
 * @brief Group.
 * @see Tides
 *
 * @author Torsten Mayer-Guerr
 * @date 2023-11-13
 *
 */
/***********************************************/

#ifndef __GROOPS_TIDESGROUP__
#define __GROOPS_TIDESGROUP__

// Latex documentation
#ifdef DOCSTRING_Tides
static const char *docstringTidesGroup = R"(
\subsection{Group}\label{tidesType:group}
Groups a set of \configClass{tides}{tidesType} and has no further effect itself.
)";
#endif

/***********************************************/

#include "classes/tides/tides.h"

/***** CLASS ***********************************/

/** @brief Group.
* @ingroup tidesGroup
* @see Tides */
class TidesGroup : public TidesBase
{
  TidesPtr tides;
  Double   factor;

public:
  TidesGroup(Config &config);

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
