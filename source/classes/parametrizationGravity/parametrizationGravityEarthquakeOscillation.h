/***********************************************/
/**
* @file parametrizationGravityEarthquakeOscillation.h
*
* @brief Earthquake oscillation.
* @see ParametrizationGravity
*
* @author Saniya Behzadpour
* @date 2017-05-08
*
*/
/***********************************************/

#ifndef __GROOPS_PARAMETRIZATIONGRAVITYEARTHQUAKE__
#define __GROOPS_PARAMETRIZATIONGRAVITYEARTHQUAKE__

// Latex documentation
#ifdef DOCSTRING_ParametrizationGravity
static const char *docstringParametrizationGravityEarthquakeOscillation = R"(
\subsection{EarthquakeOscillation}
This class is used to estimate the earthquake oscillation function parameters,
i.e. $C_{nlm}$, $\omega_{nlm}$, and $P_{nlm}$.
The results describes the variation in the gravitational potential field caused by large earthquakes.
\begin{equation}
C_{lm}(\M t) = \sum_{n=0}^NC_{nlm}(1-\cos(\omega_{nlm}d\M t)\exp(P_{nlm}\omega_{nlm}d\M t)),
\end{equation}
with $\omega_{nlm}=\frac{2\pi}{T_{nlm}}$ and $P_{nlm}=\frac{-1}{2Q_{nlm}}$ . In this equation, $Q_{nlm}$ is the attenuation factor,
$n$ is the overtone factor, $m$ is degree, $l$ is order, and $t$ is time after earthquake in second.
)";
#endif

/***********************************************/

#include "classes/sphericalHarmonicsNumbering/sphericalHarmonicsNumbering.h"
#include "classes/parametrizationGravity/parametrizationGravity.h"

/***** CLASS ***********************************/

/** @brief Earthquake oscillation.
* @ingroup parametrizationGravityGroup
* @see ParametrizationGravity */
class ParametrizationGravityEarthquakeOscillation : public ParametrizationGravityBase
{
  SphericalHarmonicsNumberingPtr numbering;
  Matrix mx;

  std::vector<std::vector<UInt>> idxC, idxS;
  Time     time0;
  UInt     _parameterCount;
  UInt     maxDegree, minDegree;
  Double   GM,R;

public:
  ParametrizationGravityEarthquakeOscillation(Config &config);

  UInt parameterCount() const {return _parameterCount;}
  void coefficients(const Time &time,  MatrixSliceRef B,  MatrixSliceRef A) const;
  void parameterName(std::vector<ParameterName> &name) const;
  void field          (const Time &time, const Vector3d &point, const Kernel &kernel, MatrixSliceRef A) const;
  void potential      (const Time &time, const Vector3d &point, MatrixSliceRef A) const;
  void radialGradient (const Time &time, const Vector3d &point, MatrixSliceRef A) const;
  void gravity        (const Time &time, const Vector3d &point, MatrixSliceRef A) const;
  void gravityGradient(const Time &time, const Vector3d &point, MatrixSliceRef A) const;
  void deformation    (const Time &time, const Vector3d &point, Double gravity, const Vector &hn, const Vector &ln, MatrixSliceRef A) const;
  SphericalHarmonics sphericalHarmonics(const Time &time, const Vector &x, UInt maxDegree) const;
  SphericalHarmonics sphericalHarmonics(const Time &time, const Vector &x, const Vector &sigma2x, UInt maxDegree) const;
};

/***********************************************/

#endif
