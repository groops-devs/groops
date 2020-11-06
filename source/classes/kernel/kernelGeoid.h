/***********************************************/
/**
* @file kernelGeoid.h
*
* @brief Integral kernel of geoid computation.
* (= Poisson Kern * gamma).
* @see Kernel
*
* @author Torsten Mayer-Guerr
* @date 2003-09-20
*
*/
/***********************************************/

#ifndef __GROOPS_KERNELGEOID__
#define __GROOPS_KERNELGEOID__

// Latex documentation
#ifdef DOCSTRING_Kernel
static const char *docstringKernelGeoid = R"(
\subsection{GeoidHeight}\label{kernelType:geoidHeight}
The geoid height is defined by Bruns formula
\begin{equation}
N = \frac{1}{\gamma}T
\end{equation}
with $T$ the disturbance potential and the normal gravity
\begin{equation}\label{normalgravity}
\gamma  = \gamma_0 - 0.30877\cdot 10^{-5}/s^2(1-0.00142\sin^2(B))h
\end{equation}
and
\begin{equation}
\gamma_0 = 9.780327\,m/s^2(1+0.0053024\sin^2(B)-0.0000058\sin^2(2B))
\end{equation}
where $h$ is the ellipsoidal height in meter and $B$ the longitude.

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

/** @brief Integral kernel of geoid computation.
* @ingroup kernelGroup
* (= Poisson Kern * gamma).
* @see Kernel */
class KernelGeoid : public Kernel
{
public:
  KernelGeoid(Config &/*config*/) {}
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
