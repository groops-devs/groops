/***********************************************/
/**
* @file kernelWaterHeight.h
*
* @brief Kernel for computation of equivalent water height.
* Conversion of density in equivalent water heights.
* @see Kernel
*
* @author Torsten Mayer-Guerr
* @date 2003-09-20
*
*/
/***********************************************/

#ifndef __GROOPS_KERNELWATERHEIGHT__
#define __GROOPS_KERNELWATERHEIGHT__

// Latex documentation
#ifdef DOCSTRING_Kernel
static const char *docstringKernelWaterHeight = R"(
\subsection{WaterHeight}\label{kernelType:waterHeight}
Height of equivalent water columns taking the effect of the loading into account.

The coefficients of the kernel defined in \eqref{eq.kernel} are
\begin{equation}
k_n = 4\pi G \rho R\frac{1+k_n'}{2n+1},
\end{equation}
where $G$ is the gravitational constant, $\rho$ is the \config{density} of water
and $k_n'$ are the load Love numbers.
)";
#endif

/***********************************************/

#include "classes/kernel/kernel.h"

/***** CLASS ***********************************/

/** @brief Kernel for computation of equivalent water height.
* @ingroup kernelGroup
* Conversion of density in equivalent water heights.
* @see Kernel */
class KernelWaterHeight : public Kernel
{
  Double rho;
  Vector kn;

public:
  KernelWaterHeight(Config &config);

  Vector coefficients       (const Vector3d &p, UInt degree) const;
  Vector inverseCoefficients(const Vector3d &p, UInt degree, Bool interior) const;
  Double inverseKernel(const Vector3d &p, const Vector3d &q, const Kernel &kernel) const;
  Double inverseKernel(const Time &time, const Vector3d &p, const GravityfieldBase &field) const;
};

/***********************************************/

#endif
