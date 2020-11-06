/***********************************************/
/**
* @file gravityfieldInInterval.h
*
* @brief Gravity fields which is not zero during a specific time span only.
* @see Gravityfield
*
* @author Torsten Mayer-Guerr
* @date 2017-11-06
*
*/
/***********************************************/

#ifndef __GROOPS_GRAVITYFIELDININTERVAL__
#define __GROOPS_GRAVITYFIELDININTERVAL__

// Latex documentation
#ifdef DOCSTRING_Gravityfield
static const char *docstringGravityfieldInInterval = R"(
\subsection{InInterval}
A \configClass{gravityfield}{gravityfieldType} is only evaluated in the interval between
\config{timeStart} inclusively and \config{timeEnd} exclusively.
Outside the interval the result is zero.

This class is useful to get a time series of monthly mean GRACE gravity field solutions.
In each month another file of potentialCoefficients is valid.
This can easily be created with \configClass{loop}{loopType}.
)";
#endif

/***********************************************/

#include "classes/gravityfield/gravityfield.h"

/***** CLASS ***********************************/

/** @brief Gravity fields which is not zero during a specific time span only.
* @ingroup gravityfieldGroup
* @see Gravityfield */
class GravityfieldInInterval : public GravityfieldBase
{
  GravityfieldPtr gravityfield;
  Time            timeStart, timeEnd;

public:
  GravityfieldInInterval(Config &config);

  Double   potential      (const Time &time, const Vector3d &point) const;
  Double   radialGradient (const Time &time, const Vector3d &point) const;
  Double   field          (const Time &time, const Vector3d &point, const Kernel &kernel) const;
  Vector3d gravity        (const Time &time, const Vector3d &point) const;
  Tensor3d gravityGradient(const Time &time, const Vector3d &point) const;
  Vector3d deformation    (const Time &time, const Vector3d &point, Double gravity, const Vector &hn, const Vector &ln) const;
  void     deformation    (const std::vector<Time> &time, const std::vector<Vector3d> &point, const std::vector<Double> &gravity,
                           const Vector &hn, const Vector &ln, std::vector< std::vector<Vector3d> > &disp) const;
  SphericalHarmonics sphericalHarmonics(const Time &time, UInt maxDegree, UInt minDegree, Double GM, Double R) const;
  void variance(const Time &time, const std::vector<Vector3d> &point, const Kernel &kernel, Matrix &D) const;
};

/***********************************************/

#endif
