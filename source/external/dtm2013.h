/***********************************************/
/**
* @file dtm2013.h
*
* @brief Fortran Wrapper.
*
* @author Torsten Mayer-Guerr
* @date   2020-02-19
*
*/
/***********************************************/

#ifndef __GROOPS_DTM2013__
#define __GROOPS_DTM2013__

#include "external/fortran.h"

struct dtm_date
{
  F77Int    type_flag;
  F77Double mjd2000;
  F77Int    day;
  F77Int    month;
  F77Int    year;
  F77Int    hour;
  F77Int    minute;
  F77Double second;
};

extern "C"
{
  void load_config(const char filePath[200]);
  void dtm_wrapper(struct dtm_date &in_date, const F77Float &alti, const F77Float &alat, const F77Float &xlon, F77Float &tz, F77Float &tinf, F77Float &tp120, F77Float &ro, F77Float &ro_unc, F77Float d[6], F77Float &wmm);
}

/***********************************************/

#endif
