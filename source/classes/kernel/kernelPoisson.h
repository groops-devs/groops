/***********************************************/
/**
* @file kernelPoisson.h
*
* @brief Poisson kernel (Potential).
* @see Kernel
*
* @author Torsten Mayer-Guerr
* @date 2003-09-20
*
*/
/***********************************************/

#ifndef __GROOPS_KERNELPOISSON__
#define __GROOPS_KERNELPOISSON__

// Latex documentation
#ifdef DOCSTRING_Kernel
static const char *docstringKernelPoisson = R"(
\subsection{Potential}
The Abel-Poisson kernel is given by
\begin{equation}
K(\cos\psi,r,R) = \frac{R(r^2-R^2)}{l^3},
\end{equation}
and the coefficients in \eqref{eq.kernel} are
\begin{equation}
k_n = 1.
\end{equation}
)";
#endif

/***********************************************/

#include "classes/kernel/kernel.h"

/***** CLASS ***********************************/

/** @brief Poisson kernel (Potential).
* @ingroup kernelGroup
* @see Kernel */
class KernelPoisson : public Kernel
{
public:
  KernelPoisson(Config &/*config*/) {}
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
