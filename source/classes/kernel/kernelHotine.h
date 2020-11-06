/***********************************************/
/**
* @file kernelHotine.h
*
* @brief Hotine kernel (gravity disturbances).
* @see Kernel
*
* @author Torsten Mayer-Guerr
* @date 2003-09-20
*
*/
/***********************************************/

#ifndef __GROOPS_KERNELHOTINE__
#define __GROOPS_KERNELHOTINE__

// Latex documentation
#ifdef DOCSTRING_Kernel
static const char *docstringKernelHotine = R"(
\subsection{Disturbance}\label{kernelType:disturbance}
Gravity disturbances in linearized form are defined by
\begin{equation}
\delta g = -\frac{dT}{dr}.
\end{equation}
The Hotine kernel is given by
\begin{equation}
K(\cos\psi,r,R) = \frac{2R^2}{l}-R\ln\frac{l+R-r\cos\psi}{r(1-\cos\psi)},
\end{equation}
and the coefficients in \eqref{eq.kernel} are
\begin{equation}
k_n = \frac{R}{n+1}.
\end{equation}
)";
#endif

/***********************************************/

#include "classes/kernel/kernel.h"

/***** CLASS ***********************************/

/** @brief Hotine kernel (gravity disturbances).
* @ingroup kernelGroup
* @see Kernel */
class KernelHotine : public Kernel
{
public:
  KernelHotine(Config &/*config*/) {}
  Double   kernel             (const Vector3d &p, const Vector3d &q) const;
  Double   radialDerivative   (const Vector3d &p, const Vector3d &q) const;
  Double   inverseKernel      (const Vector3d &p, const Vector3d &q, const Kernel &kernel) const;
  Double   inverseKernel      (const Time &time, const Vector3d &p, const GravityfieldBase &field) const;
  Vector   coefficients       (const Vector3d &p, UInt degree) const;
  Vector   inverseCoefficients(const Vector3d &p, UInt degree, Bool interior) const;
};

/***********************************************/

#endif
