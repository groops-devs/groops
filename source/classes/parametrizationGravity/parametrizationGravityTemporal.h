/***********************************************/
/**
* @file parametrizationGravityTemporal.h
*
* @brief Time variable gravity field.
* @see ParametrizationGravity
*
* @author Torsten Mayer-Guerr
* @date 2011-12-28
*
*/
/***********************************************/

#ifndef __GROOPS_PARAMETRIZATIONGRAVITYTEMPORAL__
#define __GROOPS_PARAMETRIZATIONGRAVITYTEMPORAL__

// Latex documentation
#ifdef DOCSTRING_ParametrizationGravity
static const char *docstringParametrizationGravityTemporal = R"(
\subsection{Temporal}
The time variable potential is given by
\begin{equation}
  V(\M x,t) = \sum_i V_i(\M x)\Psi_i(t),
\end{equation}
wehre $V_i(\M x)$ is the spatial parametrization of the gravity field
and can be choosen with \configClass{parametrizationGravity}{parametrizationGravityType}.
The parametrization in time domain $\Psi_i(t)$ is selected by
\configClass{parametrizationTemporal}{parametrizationTemporalType}.
The total parameter count is the parameter count of \configClass{parametrizationTemporal}{parametrizationTemporalType}
times the parameter count of \configClass{parametrizationGravity}{parametrizationGravityType}.
)";
#endif

/***********************************************/

#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "classes/parametrizationGravity/parametrizationGravity.h"

/***** CLASS ***********************************/

/** @brief Time variable gravity field.
* @ingroup parametrizationGravityGroup
* @see ParametrizationGravity */
class ParametrizationGravityTemporal : public ParametrizationGravityBase
{
  ParametrizationGravityPtr  spatial;
  ParametrizationTemporalPtr temporal;

public:
  ParametrizationGravityTemporal(Config &config);

  Bool setInterval(const Time &timeStart, const Time &timeEnd) override;
  UInt parameterCount() const override {return spatial->parameterCount() * temporal->parameterCount();}
  void parameterName(std::vector<ParameterName> &name) const override;
  void field          (const Time &time, const Vector3d &point, const Kernel &kernel, MatrixSliceRef A) const override;
  void potential      (const Time &time, const Vector3d &point, MatrixSliceRef A) const override;
  void radialGradient (const Time &time, const Vector3d &point, MatrixSliceRef A) const override;
  void gravity        (const Time &time, const Vector3d &point, MatrixSliceRef A) const override;
  void gravityGradient(const Time &time, const Vector3d &point, MatrixSliceRef A) const override;
  void deformation    (const Time &time, const Vector3d &point, Double gravity, const Vector &hn, const Vector &ln, MatrixSliceRef A) const override;
  SphericalHarmonics sphericalHarmonics(const Time &time, const Vector &x, UInt maxDegree) const override;
  SphericalHarmonics sphericalHarmonics(const Time &time, const Vector &x, const Vector &sigma2x, UInt maxDegree) const override;
};

/***********************************************/

#endif
