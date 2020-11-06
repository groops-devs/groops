/***********************************************/
/**
* @file gravityfieldEarthquakeOscillation.h
*
* @brief Earthquake oscillation.
* @see Gravityfield
*
* @author Saniya Behzadpour
* @date 2017-05-08
*
*/
/***********************************************/

#ifndef __GROOPS_GRAVITYFIELDEARTHQUAKEOSCILLATION__
#define __GROOPS_GRAVITYFIELDEARTHQUAKEOSCILLATION__

// Latex documentation
#ifdef DOCSTRING_Gravityfield
static const char *docstringGravityfieldEarthquakeOscillation = R"(
\subsection{EarthquakeOscillation}
The given \configClass{gravityfield}{gravityfieldType} is interpreted as an oscillation function
in the gravitational potential field, caused by large earthquakes.
The result is computed at time $t$ as follows:
\begin{equation}
C_{lm}(\M t) = \sum_{n=0}^NC_{nlm}(1-\cos(\omega)\exp(\frac{-\omega}{2Q_{nlm}})),
\end{equation}
with $\omega=\frac{2\pi}{T_{nlm}}(t-t_0)$. In this equation, $Q_{nlm}$ is the attenuation factor,
$n$ is the overtone factor, $m$ is degree, $l$ is order, and $t$ is time in second.
$T_{nlm}$ and $Q_{nlm}$ are computed with the elastic Earth model or observed from the long
period record of superconducting gravimeter measurements after the earthquakes.
)";
#endif

/***********************************************/

#include "classes/gravityfield/gravityfield.h"

/***** CLASS ***********************************/

/** @brief Gravityfield as oscillation.
* @ingroup gravityfieldGroup
* @see Gravityfield */
class GravityfieldEarthquakeOscillation : public GravityfieldBase
{
  Matrix   mx;
  Time     time0;
  UInt     minDegree, maxDegree;
  Double   GM, R;
  SphericalHarmonics timeVariableCoefficients(const Time &time) const;

public:
  GravityfieldEarthquakeOscillation(Config &config);

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
