/***********************************************/
/**
* @file parametrizationGravityRadialBasis.h
*
* @brief Radial basis functions.
* @see ParametrizationGravity
*
* @author Torsten Mayer-Guerr
* @date 2003-03-10
*
*/
/***********************************************/

#ifndef __GROOPS_PARAMETRIZATIONGRAVITYRADIALBASIS__
#define __GROOPS_PARAMETRIZATIONGRAVITYRADIALBASIS__

// Latex documentation
#ifdef DOCSTRING_ParametrizationGravity
static const char *docstringParametrizationGravityRadialBasis = R"(
\subsection{RadialBasis}\label{parametrizationGravityType:radialBasis}
The potential~$V$ is represented by a sum of space localizing basis functions
\begin{equation}
  V(\M x) = \sum_i a_i \Phi(\M x, \M x_i)
\end{equation}
where $a_i$ the coefficients which has to be estimated and $\Phi$ are the basis
functions given by isotropic radial \configClass{kernel}{kernelType} functions
\begin{equation}
  \Phi(\cos\psi,r,R) = \sum_n \left(\frac{R}{r}\right)^{n+1} k_n\sqrt{2n+1}\bar{P}_n(\cos\psi).
\end{equation}
The basis functions are located on a grid~$\M x_i$ given by \configClass{grid}{gridType}.
This class can also be used to estimate point masses if \configClass{kernel}{kernelType} is set to density.

The \file{parameter names}{parameterName} are \verb|*:radialBasis.<index>.<total count>:*:*|.
)";
#endif

/***********************************************/

#include "classes/parametrizationGravity/parametrizationGravity.h"

/***** CLASS ***********************************/

/** @brief Radial basis functions.
* @ingroup parametrizationGravityGroup
* @see ParametrizationGravity */
class ParametrizationGravityRadialBasis : public ParametrizationGravityBase
{
  KernelPtr kernel;                  // basis functions
  std::vector<Vector3d> sourcePoint; // center of basis functions

public:
  ParametrizationGravityRadialBasis(Config &config);

  UInt parameterCount() const override {return sourcePoint.size();}
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
