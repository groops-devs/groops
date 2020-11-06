/***********************************************/
/**
* @file blas.h
*
* @brief BLAS Wrapper fur C.
*
* @author Torsten Mayer-Guerr
* @date 2001-12-01
*
*/
/***********************************************/

#ifndef __GROOPS_BLAS__
#define __GROOPS_BLAS__

#include "external/fortran.h"

/***** DEFINES *********************************/

#define blas_dswap   FORTRANCALL(dswap     , DSWAP     )
#define blas_dscal   FORTRANCALL(dscal     , DSCAL     )
#define blas_dcopy   FORTRANCALL(dcopy     , DCOPY     )
#define blas_daxpy   FORTRANCALL(daxpy     , DAXPY     )
#define blas_ddot    FORTRANCALL(ddot      , DDOT      )
#define blas_dgemv   FORTRANCALL(wrapdgemv , WRAPDGEMV )
#define blas_dsymv   FORTRANCALL(wrapdsymv , WRAPDSYMV )
#define blas_dtrmv   FORTRANCALL(wrapdtrmv , WRAPDTRMV )
#define blas_dtrsv   FORTRANCALL(wrapdtrsv , WRAPDTRSV )
#define blas_dgemm   FORTRANCALL(wrapdgemm , WRAPDGEMM )
#define blas_dsymm   FORTRANCALL(wrapdsymm , WRAPDSYMM )
#define blas_dtrmm   FORTRANCALL(wrapdtrmm , WRAPDTRMM )
#define blas_dtrsm   FORTRANCALL(wrapdtrsm , WRAPDTRSM )
#define blas_dsyrk   FORTRANCALL(wrapdsyrk , WRAPDSYRK )
#define blas_dsyr2k  FORTRANCALL(wrapdsyr2k, WRAPDSYR2K)

extern "C"
{
void   blas_dswap (const F77Int &n,                         const F77Double x[], const F77Int &incx, F77Double y[], const F77Int &inxy);
void   blas_dscal (const F77Int &n, const F77Double &alpha,       F77Double x[], const F77Int &incx);
void   blas_dcopy (const F77Int &n,                         const F77Double x[], const F77Int &incx, F77Double y[], const F77Int &inxy);
void   blas_daxpy (const F77Int &n, const F77Double &alpha, const F77Double x[], const F77Int &incx, F77Double y[], const F77Int &inxy);
Double blas_ddot  (const F77Int &n, const F77Double x[], const F77Int &incx, const F77Double y[], const F77Int &incy);
void   blas_dgemv (const F77Bool &trans,  const F77Int  &m, const F77Int &n, const F77Double &alpha, const F77Double A[], const F77Int &ldA, const F77Double x[], const F77Int &incx, const F77Double &beta, F77Double y[], const F77Int &incy);
void   blas_dsymv (const F77Bool &upper,  const F77Int  &n, const F77Double &alpha, const F77Double A[], const F77Int &ldA, const F77Double x[], const F77Int &incx, const F77Double &beta, F77Double y[], const F77Int &incy);
void   blas_dtrmv (const F77Bool &upper,  const F77Bool &trans,  const F77Bool &unitDiag, const F77Int &n, const F77Double A[], const F77Int &ldA, F77Double x[], const F77Int &incx);
void   blas_dtrsv (const F77Bool &upper,  const F77Bool &trans,  const F77Bool &unitDiag, const F77Int &n, const F77Double A[], const F77Int &ldA, F77Double x[], const F77Int &incx);
void   blas_dgemm (const F77Bool &transA, const F77Bool &transB, const F77Int  &m, const F77Int &n, const F77Int &k, const F77Double &alpha, const F77Double A[], const F77Int &ldA, const F77Double B[], const F77Int &ldB, const F77Double &beta, F77Double C[], const F77Int &ldC);
void   blas_dsymm (const F77Bool &left,   const F77Bool &upper,  const F77Int  &m, const F77Int &n,                   const F77Double &alpha, const F77Double A[], const F77Int &ldA, const F77Double B[], const F77Int &ldB, const F77Double &beta, F77Double C[], const F77Int &ldC);
void   blas_dtrmm (const F77Bool &left,   const F77Bool &upper,  const F77Bool &trans, const F77Bool &unitDiag, const F77Int &m, const F77Int &n, const F77Double &alpha, const F77Double A[], const F77Int &ldA, F77Double B[], const F77Int &ldB);
void   blas_dtrsm (const F77Bool &left,   const F77Bool &upper,  const F77Bool &trans, const F77Bool &unitDiag, const F77Int &m, const F77Int &n, const F77Double &alpha, const F77Double A[], const F77Int &ldA, F77Double B[], const F77Int &ldB);
void   blas_dsyrk (const F77Bool &upper,  const F77Bool &trans,  const F77Int  &n, const F77Int &k, const F77Double &alpha, const F77Double A[], const F77Int &ldA, const F77Double &beta, F77Double C[], const F77Int &ldC);
void   blas_dsyr2k(const F77Bool &upper,  const F77Bool &trans,  const F77Int  &n, const F77Int &k, const F77Double &alpha, const F77Double A[], const F77Int &ldA, const F77Double B[], const F77Int &ldB, const F77Double &beta, F77Double C[], const F77Int &ldC);
}

/***********************************************/

#endif
