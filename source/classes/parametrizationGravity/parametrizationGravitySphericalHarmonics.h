/***********************************************/
/**
* @file parametrizationGravitySphericalHarmonics.h
*
* @brief Spherical Harmonics.
* @see ParametrizationGravity
*
* @author Torsten Mayer-Guerr
* @date 2003-03-10
*
*/
/***********************************************/

#ifndef __GROOPS_PARAMETRIZATIONGRAVITYHARMONICS__
#define __GROOPS_PARAMETRIZATIONGRAVITYHARMONICS__

// Latex documentation
#ifdef DOCSTRING_ParametrizationGravity
static const char *docstringParametrizationGravitySphericalHarmonics = R"(
\subsection{SphericalHarmonics}\label{parametrizationGravityType:sphericalHarmonics}
The potential~$V$ is parametrized by a expansion of (fully normalized) spherical harmonics
\begin{equation}
V(\lambda,\vartheta,r) = \frac{GM}{R}\sum_{n=0}^\infty \sum_{m=0}^n \left(\frac{R}{r}\right)^{n+1}
  \left(c_{nm} C_{nm}(\lambda,\vartheta) + s_{nm} S_{nm}(\lambda,\vartheta)\right).
\end{equation}
You can set the range of degree~$n$ with \config{minDegree} and \config{maxDegree}.
The sorting sequence of the potential coefficients in the parameter vector can be defined by
\configClass{numbering}{sphericalHarmonicsNumberingType}.

The total count of parameters is $(n_{max}+1)^2-n_{min}^2$ and
the \file{parameter names}{parameterName} are
\begin{itemize}
\item \verb|*:sphericalHarmonics.c_<degree>_<order>:*:*|,
\item \verb|*:sphericalHarmonics.s_<degree>_<order>:*:*|.
\end{itemize}
)";
#endif

/***********************************************/

#include "classes/sphericalHarmonicsNumbering/sphericalHarmonicsNumbering.h"
#include "classes/parametrizationGravity/parametrizationGravity.h"

/***** CLASS ***********************************/

/** @brief Spherical Harmonics.
* @ingroup parametrizationGravityGroup
* @see ParametrizationGravity */
class ParametrizationGravitySphericalHarmonics : public ParametrizationGravityBase
{
  SphericalHarmonicsNumberingPtr numbering;
  std::vector<std::vector<UInt>> idxC, idxS;
  UInt     _parameterCount;
  UInt     maxDegree, minDegree;
  Double   GM, R;

public:
  ParametrizationGravitySphericalHarmonics(Config &config);

  UInt parameterCount() const override {return _parameterCount;}
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
