/***********************************************/
/**
* @file gravityfieldTrend.h
*
* @brief Gravityfield as trend.
* @see Gravityfield
*
* @author Torsten Mayer-Guerr
* @date 2007-06-15
*
*/
/***********************************************/

#ifndef __GROOPS_GRAVITYFIELDTREND__
#define __GROOPS_GRAVITYFIELDTREND__

// Latex documentation
#ifdef DOCSTRING_Gravityfield
static const char *docstringGravityfieldTrend = R"(
\subsection{Trend}\label{gravityfieldType:trend}
The given \configClass{gravityfield}{gravityfieldType} is interpreted
as trend function and the result is computed at time $t$ as follows
\begin{equation}
V(\M x,t) = \frac{t-t_0}{\Delta t}V(\M x),
\end{equation}
with $t_0$ is \config{timeStart} and $\Delta t$ is \config{timeStep}.
)";
#endif

/***********************************************/

#include "classes/gravityfield/gravityfield.h"

/***** CLASS ***********************************/

/** @brief Gravityfield as trend.
* @ingroup gravityfieldGroup
* @see Gravityfield */
class GravityfieldTrend : public GravityfieldBase
{
  Time time0;
  Time timeStep;
  GravityfieldPtr gravityfield;

  Double factor(Time time) const;

public:
  GravityfieldTrend(Config &config);

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
