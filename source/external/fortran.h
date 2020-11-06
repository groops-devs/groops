/***********************************************/
/**
* @file fortran.h
*
* @brief Defines for calling fortran subroutines.
*
* @author Torsten Mayer-Guerr
* @date 2018-04-25
*
*/
/***********************************************/

#ifndef __GROOPS_FORTRAN__
#define __GROOPS_FORTRAN__

/***** DEFINES *********************************/

/// Name conventions for Fortran subroutines.
#define FORTRANCALL(lower, upper) lower ## _ // append underscore to lower-case name
// #define FOTRANCALL(lower, upper) upper    // object references in capital letters

/***** TYPES ***********************************/

/// Data types in Fortran.
//@{
typedef char   F77Char;
typedef int    F77Int;
typedef double F77Double;
typedef float  F77Float;
typedef bool   F77Bool;
//@}

/***********************************************/

#endif /* __GROOPS__ */
