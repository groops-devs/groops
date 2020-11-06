/***********************************************/
/**
* @file iers.h
*
* @brief Fortran Wrapper.
*
* @author Torsten Mayer-Guerr
* @date 2005-04-26
*
*/
/***********************************************/

#ifndef __GROOPS_IERS__
#define __GROOPS_IERS__

#include "external/fortran.h"

#define ray       FORTRANCALL(ray,       RAY)
#define ortho_eop FORTRANCALL(ortho_eop, ORTHO_EOP)
#define pmsdnut   FORTRANCALL(pmsdnut,   PMSDNUT)
#define pmsdnut2  FORTRANCALL(pmsdnut2,  PMSDNUT2)
#define utlibr    FORTRANCALL(utlibr,    UTLIBR)

extern "C" void ray(const F77Double &mjd, F77Double &corx, F77Double &cory, F77Double &cort);
extern "C" void ortho_eop(const F77Double &mjd, F77Double dxdydt[]);
extern "C" void pmsdnut(const F77Double &mjd, F77Double dxdy[]);
extern "C" void pmsdnut2(const F77Double &mjd, F77Double dxdy[]);
extern "C" void utlibr(const F77Double &mjd, F77Double &dut1, F77Double &dlod);

/***********************************************/

#endif
