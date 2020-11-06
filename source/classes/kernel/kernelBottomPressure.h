/************************************************/
/**
* @file kernelBottomPressure.h
*
* @brief Kernel for computation of ocean bottom pressure.
*
* @see Kernel
*
* @author Andreas Kvas
* @date 2017-06-20
*
*/
/***********************************************/

#ifndef __GROOPS_KERNELBOTTOMPRESSURE__
#define __GROOPS_KERNELBOTTOMPRESSURE__

// Latex documentation
#ifdef DOCSTRING_Kernel
static const char *docstringKernelBottomPressure = R"(
\subsection{BottomPressure}
Ocean bottom pressure caused by water and atmosphere masses columns taking the effect of the loading into account.

The coefficients of the kernel defined in \eqref{eq.kernel} are
\begin{equation}
k_n = \frac{4\pi G R }{\gamma}\frac{1+k_n'}{2n+1},
\end{equation}
where $G$ is the gravitational constant, $\gamma$ is the normal gravity and $k_n'$ are the load Love numbers.
)";
#endif

/***********************************************/

#include "classes/kernel/kernel.h"

/***** CLASS ***********************************/

/** @brief Kernel for computation of ocean bottom pressure.
* @ingroup kernelGroup
* Conversion of density into ocean bottom pressure.
* @see Kernel */
class KernelBottomPressure : public Kernel
{
  Vector kn;

public:
  KernelBottomPressure(Config &config);

  Vector coefficients       (const Vector3d &p, UInt degree) const;
  Vector inverseCoefficients(const Vector3d &p, UInt degree, Bool interior) const;
  Double inverseKernel(const Vector3d &p, const Vector3d &q, const Kernel &kernel) const;
  Double inverseKernel(const Time &time, const Vector3d &p, const GravityfieldBase &field) const;
};

/***********************************************/

#endif
