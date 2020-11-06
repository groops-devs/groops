/***********************************************/
/**
* @file hwm.h
*
* @brief Fortran Wrapper.
*
* @author Sandro Krauss
* @date 2019-03-21
*
*/
/***********************************************/

#ifndef __GROOPS_HWM__
#define __GROOPS_HWM__

#include "external/fortran.h"

#define hwm FORTRANCALL(hwm14, HWM14)

extern "C" void hwm(const F77Int &yd, const F77Float &sec, const F77Float &alt,
                    const F77Float &lat, const F77Float &lon, const F77Float &/*stl*/,
                    const F77Float &/*f107a*/, const F77Float &/*f107*/,
                    const F77Float ap[], F77Float outf[]);

/***********************************************/

#endif
