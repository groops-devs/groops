/***********************************************/
/**
* @file kernelFilterGauss.h
*
* @brief Kernel smoothed by an Gaussian filter.
* @see Kernel
*
* @author Torsten Mayer-Guerr
* @date 2003-09-20
*
*/
/***********************************************/

#ifndef __GROOPS_KERNELFILTERGAUSS__
#define __GROOPS_KERNELFILTERGAUSS__

// Latex documentation
#ifdef DOCSTRING_Kernel
static const char *docstringKernelFilterGauss = R"(
\subsection{FilterGauss}
Another \configClass{kernel}{kernelType} is smoothed by a gauss filter
which is defined by
\begin{equation}
F(\cos\psi) = \frac{b\cdot e^{-b(1-\cos\psi)}}{1-e^{-2b}}
\end{equation}
with $b = \frac{ln(2)}{1-\cos(r/R)}$ where $r$ is the given
smoothing \config{radius} in km and $R=6378.1366$~km is the
Earth radius.
The coefficients~$k_n$ of the \config{kernel} are multiplicated by
\begin{equation}
f_n = \frac{1}{2n+1} \int_{-1}^1 F(t)\cdot \bar{P}_n(t)\,dt.
\end{equation}
)";
#endif

/***********************************************/

#include "classes/kernel/kernel.h"

/***** CLASS ***********************************/

/** @brief Kernel smoothed by an Gaussian filter.
* @ingroup kernelGroup
* @see Kernel */
class KernelFilterGauss : public Kernel
{
  KernelPtr kernel;
  Double    radius;
  mutable Vector filterCoeff;

  Vector computeFilterCoefficients(UInt degree) const;

public:
  KernelFilterGauss(Config &config);

  Vector   coefficients       (const Vector3d &p, UInt degree) const;
  Vector   inverseCoefficients(const Vector3d &p, UInt degree, Bool interior) const;
};

/***********************************************/

#endif
