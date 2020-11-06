/***********************************************/
/**
* @file kernelSingleLayer.h
*
* @brief Single layer density.
* @see Kernel
*
* @author Torsten Mayer-Guerr
* @date 2003-09-20
*
*/
/***********************************************/

#ifndef __GROOPS_KERNELSINGLELAYER__
#define __GROOPS_KERNELSINGLELAYER__

// Latex documentation
#ifdef DOCSTRING_Kernel
static const char *docstringKernelSingleLayer = R"(
\subsection{Density}
This kernel defines a point mass or mass on a single layer ($1/l$-kernel)
taking the effect of the loading into account.

The coefficients of the kernel defined in \eqref{eq.kernel} are
\begin{equation}
k_n = 4\pi G R\frac{1+k_n'}{2n+1},
\end{equation}
where $G$ is the gravitational constant and $k_n'$ are the load Love numbers.
)";
#endif

/***********************************************/

#include "classes/kernel/kernel.h"

/***** CLASS ***********************************/

/** @brief Single layer density.
* @ingroup kernelGroup
* @see Kernel */
class KernelSingleLayer : public Kernel
{
  Vector kn;

public:
  KernelSingleLayer(Config &config);

  Double   kernel             (const Vector3d &p, const Vector3d &q) const;
  Double   radialDerivative   (const Vector3d &p, const Vector3d &q) const;
  Vector3d gradient           (const Vector3d &p, const Vector3d &q) const;
  Double   inverseKernel      (const Vector3d &p, const Vector3d &q, const Kernel &kernel) const;
  Double   inverseKernel      (const Time &time, const Vector3d &p, const GravityfieldBase &field) const;
  Vector   coefficients       (const Vector3d &p, UInt degree) const;
  Vector   inverseCoefficients(const Vector3d &p, UInt degree, Bool interior) const;
};

/***********************************************/

#endif
