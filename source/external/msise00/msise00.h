/***********************************************/
/** @file msise00.h @brief C interface to NRLMSISE-00. */
/***********************************************/
#ifndef __GROOPS_MSISE00MODEL__
#define __GROOPS_MSISE00MODEL__

#include "external/fortran.h"

extern "C" void msise00CalcWrapper(const F77Int &doy, const F77Double &sec,
                                     const F77Double &alt, const F77Double &lat,
                                     const F77Double &lon, const F77Double &lst,
                                     const F77Double &f107A, const F77Double &f107,
                                     const F77Double ap[7], F77Double &density,
                                     F77Double &temperature, F77Double &exosphericTemperature);

#endif
