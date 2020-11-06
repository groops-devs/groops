/***********************************************/
/**
* @file kernelTruncation.h
*
* @brief Kernel truncated at spherical harmonic degree
* @see Kernel
*
* @author Andreas Kvas
* @date 2019-03-01
*
*/
/***********************************************/

#ifndef __GROOPS_KERNELTRUNCATION__
#define __GROOPS_KERNELTRUNCATION__

// Latex documentation
#ifdef DOCSTRING_Kernel
static const char *docstringKernelTruncation= R"(
\subsection{Truncation}
Another \configClass{kernel}{kernelType} is truncated before \config{minDegree} and after \config{maxDegree}.
The coefficients of this kernel are defined as
\begin{equation}
  k_n =
  \begin{cases}
  1 & \text{for } n_{\text{minDegree}} \leq n \leq n_{\text{maxDegree}}\\
  0 & \text{else.} \\
  \end{cases}
\end{equation}
)";
#endif

/***********************************************/

#include "classes/kernel/kernel.h"

/***** CLASS ***********************************/

/** @brief Kernel truncated at spherical harmonic degree.
* @ingroup kernelGroup
* @see Kernel */
class KernelTruncation : public Kernel
{
  KernelPtr _kernel;
  UInt      _minDegree;
  UInt      _maxDegree;

public:
  KernelTruncation(Config &config);

  Vector   coefficients       (const Vector3d &p, UInt degree) const;
  Vector   inverseCoefficients(const Vector3d &p, UInt degree, Bool interior) const;

  UInt maxDegree() const {return _maxDegree;}
};

/***********************************************/

#endif
