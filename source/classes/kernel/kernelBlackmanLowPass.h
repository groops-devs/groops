/***********************************************/
/**
* @file kernelBlackmanLowPass.h
*
* @brief Kernel smoothed by Blackman low-pass filter
* @see Kernel
*
* @author Andreas Kvas
* @date 2019-03-01
*
*/
/***********************************************/

#ifndef __GROOPS_KERNELBLACKMANLOWPASS__
#define __GROOPS_KERNELBLACKMANLOWPASS__

// Latex documentation
#ifdef DOCSTRING_Kernel
static const char *docstringKernelBlackmanLowPass = R"(
\subsection{BlackmanLowpass}
Another \configClass{kernel}{kernelType} is smoothed by a Blackman low-pass filter. The filter is
defined through the beginning and end of the transition from pass-band to stop-band. This
transition band is specified by \config{startDegreeTransition} ($n_1$) and \config{stopDegreeTransition} ($n_2$).

The coefficients of this kernel are defined as
\begin{equation}
\begin{cases}
1 & \text{for } n < n_1 \\
A_n^2 & \text{for } n_1\leq n \leq n_2 \\
0 & \text{for } n > n_2 \\
\end{cases}
\end{equation}
with
\begin{equation}
A_n = 0.42 + 0.5\cos(\pi \frac{n-n_1}{n_2-n_1}) + 0.08 \cos(2\pi\frac{n-n_1}{n_2-n_1}).
\end{equation}
)";
#endif

/***********************************************/

#include "classes/kernel/kernel.h"

/***** CLASS ***********************************/

/** @brief Kernel smoothed by an Gaussian filter.
* @ingroup kernelGroup
* @see Kernel */
class KernelBlackmanLowPass : public Kernel
{
  KernelPtr kernel;
  UInt n1, n2;
  mutable Vector filterCoeff;

  Vector computeFilterCoefficients(UInt degree) const;

public:
  KernelBlackmanLowPass(Config &config);

  Vector   coefficients       (const Vector3d &p, UInt degree) const;
  Vector   inverseCoefficients(const Vector3d &p, UInt degree, Bool interior) const;

  UInt maxDegree() const { return n2+1; }
};

/***********************************************/

#endif
