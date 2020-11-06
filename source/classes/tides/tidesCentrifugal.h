/***********************************************/
/**
* @file tidesCentrifugal.h
*
* @brief Centrifugal force.
* Actual centrifugal force computed from Earth rotation at given time.
* @see Tides
*
* @author Torsten Mayer-Guerr
* @date 2011-10-19
*
*/
/***********************************************/

#ifndef __GROOPS_TIDESCENTRIFUGAL__
#define __GROOPS_TIDESCENTRIFUGAL__

// Latex documentation
#ifdef DOCSTRING_Tides
static const char *docstringTidesCentrifugal = R"(
\subsection{Centrifugal}\label{tidesType:centrifugal}
Computes the centrifugal potential in a rotating system
\begin{equation}
V(\M r, t) = \frac{1}{2} (\M\omega(t)\times\M r)^2.
\end{equation}
The current rotation vector $\M\omega(t)$ is computed from the
\configClass{earthRotation}{earthRotationType}
provided by the calling program.
The computed result is multiplied with \config{factor}.

Be careful, the centrifugal potential is not harmonic.
Convolution with a harmonic kernel (e.g. to compute gravity
anomalies) is not meaningful.
)";
#endif

/***********************************************/

#include "classes/tides/tides.h"

/***** CLASS ***********************************/

/** @brief Centrifugal force.
* @ingroup tidesGroup
* Actual centrifugal force computed from Earth rotation at given time.
* @see Tides */
class TidesCentrifugal : public TidesBase
{
  Double factor;

public:
  TidesCentrifugal(Config &config);

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
