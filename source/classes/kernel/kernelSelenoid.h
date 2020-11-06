/***********************************************/
/**
* @file kernelSelenoid.h
*
* @brief Integral kernel of selenoid computation.
* (= Poisson Kern * gamma).
* @see Kernel
*
* @author Beate Klinger
* @date 2013-xx-xx
*
*/
/***********************************************/

#ifndef __GROOPS_KERNELSELENOID__
#define __GROOPS_KERNELSELENOID__

// Latex documentation
#ifdef DOCSTRING_Kernel
static const char *docstringKernelSelenoid = R"(
\subsection{SelenoidHeight}
The selenoid height is defined by Bruns formula
\begin{equation}
N = \frac{1}{\gamma}T
\end{equation}
with $T$ the disturbance potential and the normal gravity $\gamma=\frac{GM}{R^2}$ of the moon.

The kernel is given by
\begin{equation}
K(\cos\psi,r,R) = \gamma\frac{R(r^2-R^2)}{l^3},
\end{equation}
and the coefficients in \eqref{eq.kernel} are
\begin{equation}
k_n = \gamma.
\end{equation}
)";
#endif

/***********************************************/

#include "classes/kernel/kernel.h"

/***** CLASS ***********************************/

/** @brief Integral kernel of selenoid computation.
* @ingroup kernelGroup
* (= Poisson Kern * gamma).
* @see Kernel */
class KernelSelenoid : public Kernel
{
public:
  KernelSelenoid(Config &/*config*/) {}
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
