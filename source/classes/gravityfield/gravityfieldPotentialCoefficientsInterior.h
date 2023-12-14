/***********************************************/
/**
* @file gravityfieldPotentialCoefficientsInterior.h
*
* @brief Potential coefficients (SphericalHarmonics).
* All classes fo this type are added together before evaluation to speed up the computation.
* @see Gravityfield
*
* @author Torsten Mayer-Guerr
* @date 2013-10-17
*
*/
/***********************************************/

#ifndef __GROOPS_GRAVITYFIELDPOTENTIALCOEFFICIENTSINTERIOR__
#define __GROOPS_GRAVITYFIELDPOTENTIALCOEFFICIENTSINTERIOR__

// Latex documentation
#ifdef DOCSTRING_Gravityfield
static const char *docstringGravityfieldPotentialCoefficientsInterior = R"(
\subsection{PotentialCoefficientsInterior}
Reads coefficients of a spherical harmonics expansion (for inner space) from file.
If set the expansion is limited in the range between \config{minDegree}
and \config{maxDegree} inclusivly. The computed result is multiplied with \config{factor}.
If \config{setSigmasToZero} is true the variances are set to zero.
This option is only important for error propagation
and does not change the result of the gravity field functionals.
)";
#endif

/***********************************************/

#include "classes/gravityfield/gravityfield.h"

/***** CLASS ***********************************/

/** @brief Potential coefficients (SphericalHarmonics).
* @ingroup gravityfieldGroup
* All classes fo this type are added together before evaluation to speed up the computation.
* @see Gravityfield */
class GravityfieldPotentialCoefficientsInterior : public GravityfieldBase
{
  SphericalHarmonics harmonics;

public:
  GravityfieldPotentialCoefficientsInterior(Config &config);
  void addPotentialCoefficients(const SphericalHarmonics &harm) {harmonics += harm;}

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

#endif
