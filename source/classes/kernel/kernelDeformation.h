/***********************************************/
/**
* @file kernelDeformation.h
*
* @brief Radial deformation by loading.
* @see Kernel
*
* @author Torsten Mayer-Guerr
* @date 2009-07-29
*
*/
/***********************************************/

#ifndef __GROOPS_KERNELDEFORMATION__
#define __GROOPS_KERNELDEFORMATION__

// Latex documentation
#ifdef DOCSTRING_Kernel
static const char *docstringKernelDeformation = R"(
\subsection{Deformation}
Computes the radial deformation caused by loading.

The coefficients of the kernel defined in \eqref{eq.kernel} are
\begin{equation}
k_n = \gamma\frac{1+k_n'}{h_n'},
\end{equation}
where $\gamma$ is the normal gravity defined in \eqref{normalgravity},
$h_n'$ and $k_n'$ are the load Love numbers and the load deformation Love numbers.
)";
#endif

/***********************************************/

#include "classes/kernel/kernel.h"

/***** CLASS ***********************************/

/** @brief Radial deformation by loading.
* @ingroup kernelGroup
* @see Kernel */
class KernelDeformation : public Kernel
{
  Matrix love;

public:
  KernelDeformation(Config &config);

  Vector   coefficients       (const Vector3d &p, UInt degree) const;
  Vector   inverseCoefficients(const Vector3d &p, UInt degree, Bool interior) const;

  UInt maxDegree() const {return love.rows()-1;}
};

/***********************************************/

#endif
