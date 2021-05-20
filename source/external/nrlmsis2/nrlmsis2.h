/***********************************************/
/**
* @file nrlmsis2.h
*
* @brief Fortran Wrapper.
*
* @author Sandro Krauss
* @date   2021-03-18
*
*/
/***********************************************/

#ifndef __GROOPS_NRLMSIS2MODEL__
#define __GROOPS_NRLMSIS2MODEL__

extern "C" void msisinitWrapper(const char *path, const char *file);
extern "C" void msiscalcWrapper(const F77Float &day, const F77Float &utsec, const F77Float &z,const F77Float &lat, const F77Float &lon, const F77Float &sfluxavg, const F77Float &sflux, const F77Float ap[7], F77Float &tn, F77Float dn[10], F77Float &tex);

/***********************************************/

#endif
