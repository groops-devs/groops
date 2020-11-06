/***********************************************/
/**
* @file kernelStokes.h
*
* @brief Stokes kernel (gravity anomalies).
* @see Kernel
*
* @author Torsten Mayer-Guerr
* @date 2003-09-20
*
*/
/***********************************************/

#ifndef __GROOPS_KERNELSTOKES__
#define __GROOPS_KERNELSTOKES__

// Latex documentation
#ifdef DOCSTRING_Kernel
static const char *docstringKernelStokes = R"(
\subsection{Anomalies}
Gravity anomalies in linearized form are defined by
\begin{equation}
\Delta g = -\frac{\partial T}{\partial r}-\frac{2}{r}T.
\end{equation}
The Stokes kernel is given by
\begin{equation}
K(\cos\psi,r,R) = \frac{2R^2}{l}-3\frac{Rl}{r^2}-\frac{R^2}{r^2}\cos\psi
\left(5+3\ln\frac{l+r-R\cos\psi}{2r}\right),
\end{equation}
and the coefficients in \eqref{eq.kernel} are
\begin{equation}
k_n = \frac{R}{n-1}.
\end{equation}
)";
#endif

/***********************************************/

#include "classes/kernel/kernel.h"

/***** CLASS ***********************************/

/** @brief Stokes kernel (gravity anomalies).
* @ingroup kernelGroup
* @see Kernel */
class KernelStokes : public Kernel
{
public:
  KernelStokes(Config &/*config*/) {}
  Double   kernel             (const Vector3d &p, const Vector3d &q) const;
  Double   radialDerivative   (const Vector3d &p, const Vector3d &q) const;
  Vector3d gradient           (const Vector3d &p, const Vector3d &q) const;
  Tensor3d gradientGradient   (const Vector3d &p, const Vector3d &q) const;
  Double   inverseKernel      (const Vector3d &p, const Vector3d &q, const Kernel &kernel) const;
  Double   inverseKernel      (const Time &time, const Vector3d &p, const GravityfieldBase &field) const;
  Vector   coefficients       (const Vector3d &p, UInt degree) const;
  Vector   inverseCoefficients(const Vector3d &p, UInt degree, Bool interior) const;
};

/***********************************************/

#endif
