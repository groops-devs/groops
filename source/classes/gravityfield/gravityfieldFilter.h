/***********************************************/
/**
* @file gravityfieldFilter.h
*
* @brief Filtered spherical harmonics.
* @see Gravityfield
*
* @author Torsten Mayer-Guerr
* @date 2020-07-20
*
*/
/***********************************************/

#ifndef __GROOPS_GRAVITYFIELDFILTER__
#define __GROOPS_GRAVITYFIELDFILTER__

// Latex documentation
#ifdef DOCSTRING_Gravityfield
static const char *docstringGravityfieldFilter = R"(
\subsection{Filter}
Convert \configClass{gravityfield}{gravityfieldType} to spherical harmonics
and \configClass{filter}{sphericalHarmonicsFilterType} the coefficients.
)";
#endif

/***********************************************/

#include "classes/sphericalHarmonicsFilter/sphericalHarmonicsFilter.h"
#include "classes/gravityfield/gravityfield.h"

/***** CLASS ***********************************/

/** @brief Filtered spherical harmonics.
* @ingroup gravityfieldGroup
* @see Gravityfield */
class GravityfieldFilter : public GravityfieldBase
{
  GravityfieldPtr             gravityfield;
  SphericalHarmonicsFilterPtr filter;

public:
  GravityfieldFilter(Config &config);

  Double   potential      (const Time &time, const Vector3d &point) const;
  Double   radialGradient (const Time &time, const Vector3d &point) const;
  Double   field          (const Time &time, const Vector3d &point, const Kernel &kernel) const;
  Vector3d gravity        (const Time &time, const Vector3d &point) const;
  Tensor3d gravityGradient(const Time &time, const Vector3d &point) const;
  Vector3d deformation    (const Time &time, const Vector3d &point, Double gravity, const Vector &hn, const Vector &ln) const;
  void     deformation    (const std::vector<Time> &time, const std::vector<Vector3d> &point, const std::vector<Double> &gravity,
                           const Vector &hn, const Vector &ln, std::vector<std::vector<Vector3d>> &disp) const;

  SphericalHarmonics sphericalHarmonics(const Time &time, UInt maxDegree=INFINITYDEGREE, UInt minDegree=0, Double GM=0.0, Double R=0.0) const;

  void   variance(const Time &time, const std::vector<Vector3d> &point, const Kernel &kernel, Matrix &D) const;
  Double variance(const Time &time, const Vector3d &point, const Kernel &kernel) const;
};

/***********************************************/

#endif /* __GROOPS_GRAVITYFIELD__ */
