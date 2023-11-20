/***********************************************/
/**
* @file parametrizationGravityLinearTransformation.h
*
* @brief Gravity field parametrization based on the linear transformation of another parametrizationGravity.
* @see ParametrizationGravity
*
* @author Andreas Kvas
* @date 2020-06-28
*
*/
/***********************************************/

#ifndef __GROOPS_PARAMETRIZATIONGRAVITYLINEARTRANSFORMATION__
#define __GROOPS_PARAMETRIZATIONGRAVITYLINEARTRANSFORMATION__

// Latex documentation
#ifdef DOCSTRING_ParametrizationGravity
static const char *docstringParametrizationGravityLinearTransformation = R"(
\subsection{LinearTransformation}
Parametrization of the gravity field on the basis of a linear transformation of a source parametrization.
The linear transformation changes the original solution space represented by
\configClass{pararametrizationGravitySource}{parametrizationGravityType} from
\begin{equation}
  \mathbf{l} = \mathbf{A}\mathbf{x} + \mathbf{e}
\end{equation}
to
\begin{equation}
  \mathbf{l} = \mathbf{A}\mathbf{F}\mathbf{y} + \mathbf{e}
\end{equation}
through the linear transformation $\mathbf{x}=\mathbf{F}\mathbf{y}$.
It follows that the rows of the matrix $\mathbf{F}$ in \file{inputfileTransformationMatrix}{matrix} coincides with
the number of parameters in \configClass{pararametrizationGravitySource}{parametrizationGravityType}.
The new parameter count is given by the number of columns in $\mathbf{F}$ and may be smaller, equal or larger
than the original parameter count.

The \file{parameter names}{parameterName} are \verb|*:transformedParameter.<index>.<total count>:*:*|.
)";
#endif

/***********************************************/

#include "classes/parametrizationGravity/parametrizationGravity.h"

/***** CLASS ***********************************/

/** @brief Gravity field parametrization based on the linear transformation of another parametrizationGravity.
* @ingroup parametrizationGravityGroup
* @see ParametrizationGravity */
class ParametrizationGravityLinearTransformation : public ParametrizationGravityBase
{
  ParametrizationGravityPtr  source;
  Matrix                     transformationMatrix;

public:
  ParametrizationGravityLinearTransformation(Config &config);

  UInt parameterCount() const override {return transformationMatrix.columns();}
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
