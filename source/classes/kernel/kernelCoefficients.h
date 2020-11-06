/***********************************************/
/**
* @file kernelCoefficients.h
*
* @brief Kernel from coefficients.
* Given as series of Legendre polynomials.
* @see Kernel
*
* @author Torsten Mayer-Guerr
* @date 2003-09-20
*
*/
/***********************************************/

#ifndef __GROOPS_KERNELCOEFF__
#define __GROOPS_KERNELCOEFF__

// Latex documentation
#ifdef DOCSTRING_Kernel
static const char *docstringKernelCoefficients = R"(
\subsection{Coefficients}\label{kernelType:coefficients}
The kernel is defined by the coefficients $k_n$ given by file.
)";
#endif

/***********************************************/

#include "classes/kernel/kernel.h"

/***** CLASS ***********************************/

/** @brief Kernel from coefficients.
* Given as series of Legendre polynomials.
* @ingroup kernelGroup
* @see Kernel */
class KernelCoefficients : public Kernel
{
  Vector kn;

public:
  KernelCoefficients(Config &config);

  Vector   coefficients       (const Vector3d &p, UInt degree) const;
  Vector   inverseCoefficients(const Vector3d &p, UInt degree, Bool interior) const;

  UInt maxDegree() const {return kn.rows()-1;}
};

/***********************************************/

#endif

