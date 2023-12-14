/***********************************************/
/**
* @file gravityfieldPotentialCoefficients.h
*
* @brief Potential coefficients (SphericalHarmonics).
* All classes fo this type are added together before evaluation to speed up the computation.
* @see Gravityfield
*
* @author Torsten Mayer-Guerr
* @date 2001-08-25
*
*/
/***********************************************/

#ifndef __GROOPS_GRAVITYFIELDPOTENTIALCOEFFICIENTS__
#define __GROOPS_GRAVITYFIELDPOTENTIALCOEFFICIENTS__

// Latex documentation
#ifdef DOCSTRING_Gravityfield
static const char *docstringGravityfieldPotentialCoefficients = R"(
\subsection{PotentialCoefficients}\label{gravityfieldType:potentialCoefficients}
Reads coefficients of a spherical harmonics expansion from file.
The potential is given by
\begin{equation}
V(\lambda,\vartheta,r) = \frac{GM}{R}\sum_{n=0}^\infty \sum_{m=0}^n \left(\frac{R}{r}\right)^{n+1}
  \left(c_{nm} C_{nm}(\lambda,\vartheta) + s_{nm} S_{nm}(\lambda,\vartheta)\right).
\end{equation}
If set the expansion is limited in the range between \config{minDegree}
and \config{maxDegree} inclusivly. The computed result
is multiplied with \config{factor}. If \config{setSigmasToZero} is true
the variances are set to zero. This option is only important for variance propagation
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
class GravityfieldPotentialCoefficients : public GravityfieldBase
{
  SphericalHarmonics harmonics;

public:
  GravityfieldPotentialCoefficients(Config &config);
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
