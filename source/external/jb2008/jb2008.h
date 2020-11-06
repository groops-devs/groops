/***********************************************/
/**
* @file jb2008.h
*
* @brief Fortran Wrapper.
*
* @author Beate Klinger
* @date   2015-12-03
*
*/
/***********************************************/

#ifndef __GROOPS_JB2008MODEL__
#define __GROOPS_JB2008MODEL__

#include "external/fortran.h"

#define jb2008 FORTRANCALL(jb2008, JB2008)

extern "C" void jb2008(const F77Double &AMJD, const F77Double SUN[2], const F77Double SAT[3],
                       const F77Double &F10,  const F77Double &F10B,
                       const F77Double &S10,  const F77Double &S10B,
                       const F77Double &M10,  const F77Double &M10B,
                       const F77Double &Y10,  const F77Double &Y10B,
                       const F77Double &DSTDTC,
                       F77Double TEMP[2], F77Double &RHO);

/***********************************************/

#endif
