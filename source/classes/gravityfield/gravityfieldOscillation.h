/***********************************************/
/**
* @file gravityfieldOscillation.h
*
* @brief Gravityfield as oscillation.
* @see Gravityfield
*
* @author Torsten Mayer-Guerr
* @date 2015-06-07
*
*/
/***********************************************/

#ifndef __GROOPS_GRAVITYFIELDOSCILLATION__
#define __GROOPS_GRAVITYFIELDOSCILLATION__

// Latex documentation
#ifdef DOCSTRING_Gravityfield
static const char *docstringGravityfieldOscillation = R"(
\subsection{Oscillation}\label{gravityfieldType:oscillation}
The given \configClass{gravityfield}{gravityfieldType} is interpreted
as oscillation function and the result is computed at time $t$ as follows
\begin{equation}
V(\M x,t) = \cos(\omega)V_{cos}(\M x)+\sin(\omega)V_{sin}(\M x),
\end{equation}
with $\omega=\frac{2\pi}{T}(t-t_0)$.
)";
#endif

/***********************************************/

#include "classes/gravityfield/gravityfield.h"

/***** CLASS ***********************************/

/** @brief Gravityfield as oscillation.
* @ingroup gravityfieldGroup
* @see Gravityfield */
class GravityfieldOscillation : public GravityfieldBase
{
  Time time0;
  Time timePeriod;
  GravityfieldPtr gravityfieldCos;
  GravityfieldPtr gravityfieldSin;

public:
  GravityfieldOscillation(Config &config);

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

#endif /* __GROOPS_GRAVITYFIELD__ */
