/***********************************************/
/**
* @file kernelRadialGradient.h
*
* @brief 2nd radial derivative (gravity gradients).
* @see Kernel
*
* @author Torsten Mayer-Guerr
* @date 2003-09-20
*
*/
/***********************************************/

#ifndef __GROOPS_KERNELRADIALGRADIENT__
#define __GROOPS_KERNELRADIALGRADIENT__

// Latex documentation
#ifdef DOCSTRING_Kernel
static const char *docstringKernelRadialGradient = R"(
\subsection{RadialGradient}
This kernel defines the second radial derivative of the (disturbance) potential.
\begin{equation}
T_{rr} = \frac{\partial^2 T}{\partial r^2}.
\end{equation}
The coefficients of the kernel defined in \eqref{eq.kernel} are
\begin{equation}
k_n = \frac{r^2}{(n+1)(n+2)}.
\end{equation}
)";
#endif

/***********************************************/

#include "classes/kernel/kernel.h"

/***** CLASS ***********************************/

/** @brief 2nd radial derivative (gravity gradients).
* @ingroup kernelGroup
* @see Kernel */
class KernelRadialGradient : public Kernel
{
public:
  KernelRadialGradient(Config &/*config*/) {}
  Vector   coefficients       (const Vector3d &p, UInt degree) const;
  Vector   inverseCoefficients(const Vector3d &p, UInt degree, Bool interior) const;
};

/***********************************************/

#endif
