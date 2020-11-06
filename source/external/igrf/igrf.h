/***********************************************/
/**
* @file igrf.h
*
* @brief Fortran Wrapper.
*
* @author Sebastian Strasser
* @date 2016-03-10
*
*/
/***********************************************/

#ifndef __GROOPS_IGRF__
#define __GROOPS_IGRF__

#include "external/fortran.h"

#define igrfSynthesis FORTRANCALL(igrf13syn, IGRF13SYN)

extern "C" void igrfSynthesis(const F77Int &isv, const F77Double &date, const F77Int &itype, const F77Double &alt, const F77Double &colat, const F77Double &elong, F77Double &x, F77Double &y, F77Double &z, F77Double &f);

/***********************************************/

#endif
