/***********************************************/
/**
* @file gravityfieldGroup.h
*
* @brief Group.
* @see Gravityfield
*
* @author Torsten Mayer-Guerr
* @date 2023-11-13
*
*/
/***********************************************/

#ifndef __GROOPS_GravityfieldGroup__
#define __GROOPS_GravityfieldGroup__

// Latex documentation
#ifdef DOCSTRING_Gravityfield
static const char *docstringGravityfieldGroup = R"(
\subsection{Group}\label{gravityfieldType:group}
Groups a set of \configClass{gravityfield}{gravityfieldType} and has no further effect itself.
)";
#endif

/***********************************************/

#include "classes/gravityfield/gravityfield.h"

/***** CLASS ***********************************/

/** @brief Group.
* @ingroup gravityfieldGroup
* @see Gravityfield */
class GravityfieldGroup : public GravityfieldBase
{
  GravityfieldPtr gravityfield;
  Double          factor;

public:
  GravityfieldGroup(Config &config);

  Double   potential      (const Time &time, const Vector3d &point) const;
  Double   radialGradient (const Time &time, const Vector3d &point) const;
  Double   field          (const Time &time, const Vector3d &point, const Kernel &kernel) const;
  Vector3d gravity        (const Time &time, const Vector3d &point) const;
  Tensor3d gravityGradient(const Time &time, const Vector3d &point) const;
  Vector3d deformation    (const Time &time, const Vector3d &point, Double gravity, const Vector &hn, const Vector &ln) const;
  void     deformation    (const std::vector<Time> &time, const std::vector<Vector3d> &point, const std::vector<Double> &gravity,
                           const Vector &hn, const Vector &ln, std::vector<std::vector<Vector3d>> &disp) const;

  SphericalHarmonics sphericalHarmonics(const Time &time, UInt maxDegree, UInt minDegree, Double GM, Double R) const;

  void variance(const Time &time, const std::vector<Vector3d> &point, const Kernel &kernel, Matrix &D) const;
};

/***********************************************/

#endif /* __GROOPS_GRAVITYFIELD__ */
