/***********************************************/
/**
* @file lapack.h
*
* LAPACK Wrapper.
*
* @author Torsten Mayer-Guerr
* @date 2001-12-01
*
*/
/***********************************************/

#ifndef __GROOPS_LAPACK__
#define __GROOPS_LAPACK__

#include "external/fortran.h"

/***** DEFINES *********************************/

#define wrapdlacpy  FORTRANCALL(wrapdlacpy, WRAPDLACPY)
#define wrapdpotrf  FORTRANCALL(wrapdpotrf, WRAPDPOTRF)
#define wrapdpotri  FORTRANCALL(wrapdpotri, WRAPDPOTRI)
#define wrapdpstrf  FORTRANCALL(wrapdpstrf, WRAPDPSTRF)
#define wrapdtrtri  FORTRANCALL(wrapdtrtri, WRAPDTRTRI)
#define wrapdlauum  FORTRANCALL(wrapdlauum, WRAPDLAUUM)
#define wrapdgetrf  FORTRANCALL(wrapdgetrf, WRAPDGETRF)
#define wrapdgesv   FORTRANCALL(wrapdgesv , WRAPDGESV )
#define wrapdsgesv  FORTRANCALL(wrapdsgesv, WRAPDSGESV)
#define wrapdgetri  FORTRANCALL(wrapdgetri, WRAPDGETRI)
#define wrapdgels   FORTRANCALL(wrapdgels , WRAPDGELS )
#define wrapdgeqrf  FORTRANCALL(wrapdgeqrf, WRAPDGEQRF)
#define wrapdormqr  FORTRANCALL(wrapdormqr, WRAPDORMQR)
#define wrapdorgqr  FORTRANCALL(wrapdorgqr, WRAPDORGQR)
#define wrapdpbsv   FORTRANCALL(wrapdpbsv , WRAPDPBSV )
#define wrapdpbtrf  FORTRANCALL(wrapdpbtrf, WRAPDPBTRF)
#define wrapdgbtrf  FORTRANCALL(wrapdgbtrf, WRAPDGBTRF)
#define wrapdtbtrs  FORTRANCALL(wrapdtbtrs, WRAPDTBTRS)
#define wrapdgbtrs  FORTRANCALL(wrapdgbtrs, WRAPDGBTRS)
#define wrapdsyev   FORTRANCALL(wrapdsyev,  WRAPDSYEV )
#define wrapdgeev   FORTRANCALL(wrapdgeev , WRAPDGEEV )
#define wrapdgesvd  FORTRANCALL(wrapdgesvd, WRAPDGESVD)
#define wrapdgesdd  FORTRANCALL(wrapdgesdd, WRAPDGESDD)

/***********************************************/

// copy
void lapack_dlacpy(UInt m, UInt n, const Double A[], UInt ldA, Double B[], UInt ldB);

// Cholesky, Inverse
Int lapack_dpotrf(Bool upper, UInt n, Double A[], UInt ldA);
Int lapack_dpotri(Bool upper, UInt n, Double A[], UInt ldA);
Int lapack_dtrtri(Bool upper, UInt n, Double A[], UInt ldA);
Int lapack_dlaumm(Bool upper, UInt n, Double A[], UInt ldA);

// LU-decomposition, Solve, Inverse
Int lapack_dgetrf(UInt m, UInt n, Double A[], UInt ldA, F77Int ipiv[]);
Int lapack_dgesv (UInt n, UInt nrhs, Double A[], UInt ldA, F77Int ipiv[], Double B[], UInt ldB);
Int lapack_dsgesv(UInt n, UInt nrhs, Double A[], UInt ldA, F77Int ipiv[], Double B[], UInt ldB, Double X[], UInt ldX, Int iter);
Int lapack_dgetri(UInt n, Double A[], UInt ldA, F77Int ipiv[]);
Int lapack_dpstrf(Bool upper, UInt n, Double A[], UInt ldA, F77Int ipiv[], F77Int &rank, Double tol);

// Least-Square
Int lapack_dgels (Bool trans, UInt m, UInt n, UInt nrhs, Double A[], UInt ldA, Double B[], UInt ldB);

// QR-decomposition
Int lapack_dgeqrf(UInt m, UInt n, Double A[], UInt ldA, Double tau[]);
Int lapack_dormqr(Bool left, Bool trans, UInt m, UInt n, UInt k, const Double A[], UInt ldA, const Double tau[], Double C[], UInt ldC);
Int lapack_dorgqr(UInt m, UInt n, UInt k, Double A[], UInt ldA, const Double tau[]);

// Band matrices
Int lapack_dpbsv(Bool upper, UInt n, UInt kd, UInt nrhs, Double A[], UInt ldA, Double B[], UInt ldB);
Int lapack_dpbtrf(Bool upper, UInt n, UInt kd, Double A[], UInt ldA);
Int lapack_dgbtrf(UInt n, UInt m, UInt kl, UInt ku, Double A[], UInt ldA, F77Int ipiv[]);
Int lapack_dtbtrs(Bool upper, Bool trans, Bool unitDiag, UInt n, UInt kd, UInt nrhs, const Double A[], UInt ldA, Double B[], UInt ldB);
Int lapack_dgbtrs(Bool trans, UInt n, UInt kl, UInt ku, UInt nrhs, const Double A[], UInt ldA, F77Int ipiv[], Double B[], UInt ldB);

// eigen value decomposition
Int lapack_dsyev(Bool jobz, Bool upper, UInt n, Double A[], UInt ldA, Double W[]);
Int lapack_dgeev(Bool jobvl, Bool jobvr, UInt n, Double A[], UInt ldA, Double WR[], Double WI[], Double VL[], UInt ldVL, Double VR[], UInt ldVR);

// sigular value decomposition
Int lapack_dgesvd(Bool jobu, Bool jobvt, UInt m, UInt n, Double A[], UInt ldA, Double S[], Double U[], UInt ldU, Double VT[], UInt ldvt);
Int lapack_dgesdd(Bool jobz, UInt m, UInt n, Double A[], UInt ldA, Double S[], Double U[], UInt ldU, Double VT[], UInt ldvt);

/***********************************************/
/***********************************************/

extern "C"
{
void wrapdlacpy(const F77Int &m, const F77Int &n, const F77Double A[], const F77Int &ldA, F77Double B[], const F77Int &ldB);
void wrapdpotrf(const F77Int &upper, const F77Int &n, F77Double A[], const F77Int &ldA, F77Int &info);
void wrapdpotri(const F77Int &upper, const F77Int &n, F77Double A[], const F77Int &ldA, F77Int &info);
void wrapdpstrf(const F77Int &upper, const F77Int &n, F77Double A[], const F77Int &ldA, F77Int ipiv[], F77Int &rank, F77Double &tol, F77Double work[], F77Int &info);
void wrapdtrtri(const F77Int &upper, const F77Int &n, F77Double A[], const F77Int &ldA, F77Int &info);
void wrapdlauum(const F77Int &upper, const F77Int &n, F77Double A[], const F77Int &ldA, F77Int &info);
void wrapdgetrf(const F77Int &m, const F77Int &n,    F77Double A[], const F77Int &ldA, F77Int ipiv[], F77Int &info);
void wrapdgesv (const F77Int &n, const F77Int &nrhs, F77Double A[], const F77Int &ldA, F77Int ipiv[], F77Double B[], const F77Int &ldB, F77Int &info);
void wrapdsgesv(const F77Int &n, const F77Int &nrhs, F77Double A[], const F77Int &ldA, F77Int ipiv[], F77Double B[], const F77Int &ldB, F77Double X[], const F77Int &ldX, F77Double work[], F77Float swork[], F77Int &iter, F77Int &info);
void wrapdgetri(const F77Int &n, F77Double A[], const F77Int &ldA, F77Int ipiv[], F77Double work[], const F77Int &lwork, F77Int &info);
void wrapdgels (const F77Int &trans, const F77Int &m, const F77Int &n, const F77Int &nrhs, F77Double A[], const F77Int &ldA, F77Double B[], const F77Int &ldB, F77Double work[], const F77Int &lwork, F77Int &info);
void wrapdgeqrf(const F77Int &m, const F77Int &n, F77Double A[], const F77Int &ldA, F77Double tau[], F77Double work[], const F77Int &lwork, F77Int &info);
void wrapdormqr(const F77Int &left, const F77Int &trans, const F77Int &m, const F77Int &n, const F77Int &k, const F77Double A[], const F77Int &ldA, const F77Double tau[], F77Double C[], const F77Int &ldC, F77Double work[], const F77Int &lwork, F77Int &info);
void wrapdorgqr(const F77Int &m, const F77Int &n, const F77Int &k, F77Double A[], const F77Int &ldA, const F77Double tau[], F77Double work[], const F77Int &lwork, F77Int &info);
void wrapdpbsv (const F77Int &upper, const F77Int &n, const F77Int &kd, const F77Int &nrhs, F77Double A[], const F77Int &ldA, F77Double B[], const F77Int &ldB, F77Int &info);
void wrapdpbtrf(const F77Int &upper, const F77Int &n, const F77Int &kd, F77Double A[], const F77Int &ldA, F77Int &info);
void wrapdgbtrf(const F77Int &n, const F77Int &m, const F77Int &kl, const F77Int &ku, F77Double A[], const F77Int &ldA, F77Int ipiv[], F77Int &info);
void wrapdtbtrs(const F77Int &upper,  const F77Int &trans, const F77Int &unitDiag, const F77Int &n, const F77Int &kd, const F77Int &nrhs, const F77Double A[], const F77Int &ldA, F77Double B[], const F77Int &ldB, F77Int &info);
void wrapdgbtrs(const F77Int &trans, const F77Int &n, const F77Int &kl, const F77Int &ku, const F77Int &nrhs, const F77Double A[], const F77Int &ldA, F77Int ipiv[], F77Double B[], const F77Int &ldB, F77Int &info);
void wrapdsyev (const F77Int &jobz, const F77Int &upper, const F77Int &n, F77Double A[], const F77Int &ldA, F77Double W[], F77Double work[], const F77Int &lwork, F77Int &info);
void wrapdgeev (const F77Int &jobvl, const F77Int &jobvr, const F77Int &n, F77Double A[], const F77Int &ldA, F77Double WR[], F77Double WI[], F77Double VL[], const F77Int &ldVL, F77Double VR[], const F77Int &ldVR, F77Double work[], const F77Int &lwork, F77Int &info);
void wrapdgesvd(const F77Int &jobu, const F77Int &jobvt, const F77Int &m, const F77Int &n, F77Double A[], const F77Int &ldA, F77Double S[], F77Double U[], const F77Int &ldU, F77Double VT[], const F77Int &ldvt, F77Double work[], const F77Int &lwork, F77Int &info);
void wrapdgesdd(const F77Int &jobz, const F77Int &m, const F77Int &n, F77Double A[], const F77Int &ldA, F77Double S[], F77Double U[], const F77Int &ldU, F77Double VT[], const F77Int &ldvt, F77Double work[], const F77Int &lwork, F77Int iwork[], F77Int &info);
}

/***** INLINES *********************************/

inline void lapack_dlacpy(UInt m, UInt n, const Double A[], UInt ldA, Double B[], UInt ldB)
{
  wrapdlacpy(static_cast<F77Int>(m), static_cast<F77Int>(n), A, static_cast<F77Int>(ldA), B, static_cast<F77Int>(ldB));
}

inline Int lapack_dpotrf(Bool upper, UInt n, Double A[], UInt ldA)
{
  F77Int info;
  wrapdpotrf(upper, static_cast<F77Int>(n), A, static_cast<F77Int>(ldA), info);
  return info;
}

inline Int lapack_dpotri(Bool upper, UInt n, Double A[], UInt ldA)
{
  F77Int info;
  wrapdpotri(upper, static_cast<F77Int>(n), A,static_cast<F77Int>(ldA), info);
  return info;
}

inline Int lapack_dpstrf(Bool upper, UInt n, Double A[], UInt ldA, F77Int ipiv[], F77Int &rank, Double tol)
{
  F77Int info;
  F77Double *work  = new F77Double[2*n];
  wrapdpstrf(upper, static_cast<F77Int>(n), A, static_cast<F77Int>(ldA), ipiv, rank, tol, work,info);
  delete[] work;
  return info;
}

inline Int lapack_dtrtri(Bool upper, UInt n, Double A[], UInt ldA)
{
  F77Int info;
  wrapdtrtri(upper, static_cast<F77Int>(n), A, static_cast<F77Int>(ldA), info);
  return info;
}

inline Int lapack_dlauum(Bool upper, UInt n, Double A[], UInt ldA)
{
  F77Int info;
  wrapdlauum(upper, static_cast<F77Int>(n), A, static_cast<F77Int>(ldA), info);
  return info;
}

inline Int lapack_dgetrf(UInt m, UInt n, Double A[], UInt ldA, F77Int ipiv[])
{
  F77Int info;
  wrapdgetrf(static_cast<F77Int>(m), static_cast<F77Int>(n), A, static_cast<F77Int>(ldA), ipiv, info);
  return info;
}

inline Int lapack_dgesv (UInt n, UInt nrhs, Double A[], UInt ldA, F77Int ipiv[], Double B[], UInt ldB)
{
  F77Int info;
  wrapdgesv(static_cast<F77Int>(n), static_cast<F77Int>(nrhs), A, static_cast<F77Int>(ldA), ipiv, B, static_cast<F77Int>(ldB), info);
  return info;
}

inline Int lapack_dsgesv (UInt n, UInt nrhs, Double A[], UInt ldA, F77Int ipiv[], Double B[], UInt ldB, Double X[], UInt ldX, Int iter)
{
  F77Int info;
  F77Double *work  = new F77Double[n*nrhs];
  F77Float  *swork = new F77Float [n*(n+nrhs)];
  wrapdsgesv(static_cast<F77Int>(n), static_cast<F77Int>(nrhs), A, static_cast<F77Int>(ldA), ipiv, B, static_cast<F77Int>(ldB), X, static_cast<F77Int>(ldX), work, swork, iter, info);
  delete[] work;
  delete[] swork;
  return info;
}

inline Int lapack_dgetri(UInt n, Double A[], UInt ldA, F77Int ipiv[])
{
  F77Int info;
  F77Double tmp;
  wrapdgetri(static_cast<F77Int>(n), A, static_cast<F77Int>(ldA), ipiv, &tmp, -1, info);
  F77Int lwork = (tmp>0) ? static_cast<F77Int>(tmp) : static_cast<F77Int>(n);
  F77Double *work = new F77Double[lwork];
  wrapdgetri(static_cast<F77Int>(n), A, static_cast<F77Int>(ldA), ipiv, work, lwork, info);
  delete[] work;
  return info;
}

inline Int lapack_dgels (Bool trans, UInt m, UInt n, UInt nrhs, Double A[], UInt ldA, Double B[], UInt ldB)
{
  F77Int    info;
  F77Double tmp;
  wrapdgels(trans,static_cast<F77Int>(m), static_cast<F77Int>(n), static_cast<F77Int>(nrhs), A, static_cast<F77Int>(ldA), B, static_cast<F77Int>(ldB), &tmp, -1, info);
  F77Int lwork = (tmp>0) ? static_cast<F77Int>(tmp) : static_cast<F77Int>(std::max(n,m));
  F77Double *work = new F77Double[lwork];
  wrapdgels(trans, static_cast<F77Int>(m), static_cast<F77Int>(n), static_cast<F77Int>(nrhs), A, static_cast<F77Int>(ldA), B, static_cast<F77Int>(ldB), work, lwork, info);
  delete[] work;
  return info;
}

inline Int lapack_dgeqrf(UInt m, UInt n, Double A[], UInt ldA, Double tau[])
{
  F77Int    info;
  F77Double tmp;
  wrapdgeqrf(static_cast<F77Int>(m), static_cast<F77Int>(n), A, static_cast<F77Int>(ldA), tau, &tmp, -1, info);
  F77Int lwork = (tmp>0) ? static_cast<F77Int>(tmp) : static_cast<F77Int>(std::max(n,m));
  F77Double *work = new F77Double[lwork];
  wrapdgeqrf(static_cast<F77Int>(m), static_cast<F77Int>(n), A, static_cast<F77Int>(ldA), tau, work, lwork, info);
  delete[] work;
  return info;
}

inline Int lapack_dormqr(Bool left, Bool trans, UInt m, UInt n, UInt k, const Double A[], UInt ldA, const Double tau[], Double C[], UInt ldC)
{
  F77Int    info;
  F77Double tmp;
  wrapdormqr(left, trans, static_cast<F77Int>(m), static_cast<F77Int>(n), static_cast<F77Int>(k), A, static_cast<F77Int>(ldA), tau, C, static_cast<F77Int>(ldC), &tmp, -1, info);
  F77Int lwork = (tmp>0) ? static_cast<F77Int>(tmp) : static_cast<F77Int>(std::max(n,m));
  F77Double *work = new F77Double[lwork];
  wrapdormqr(left,trans,static_cast<F77Int>(m), static_cast<F77Int>(n), static_cast<F77Int>(k), A,static_cast<F77Int>(ldA), tau, C, static_cast<F77Int>(ldC), work, lwork, info);
  delete[] work;
  return info;
}

inline Int lapack_dorgqr(UInt m, UInt n, UInt k, Double A[], UInt ldA, const Double tau[])
{
  F77Int    info;
  F77Double tmp;
  wrapdorgqr(static_cast<F77Int>(m), static_cast<F77Int>(n), static_cast<F77Int>(k),A,static_cast<F77Int>(ldA), tau,&tmp,-1,info);
  F77Int lwork = (tmp>0) ? static_cast<F77Int>(tmp) : static_cast<F77Int>(std::max(n,m));

  F77Double *work = new F77Double[lwork];
  wrapdorgqr(static_cast<F77Int>(m), static_cast<F77Int>(n), static_cast<F77Int>(k),A,static_cast<F77Int>(ldA), tau,work,lwork,info);
  delete[] work;
  return info;
}

inline Int lapack_dpbsv(Bool upper, UInt n, UInt kd, UInt nrhs, Double A[], UInt ldA, Double B[], UInt ldB)
{
  F77Int info;
  wrapdpbsv(upper,static_cast<F77Int>(n), static_cast<F77Int>(kd), static_cast<F77Int>(nrhs), A,static_cast<F77Int>(ldA), B,static_cast<F77Int>(ldB), info);
  return info;
}

inline Int lapack_dpbtrf(Bool upper, UInt n, UInt kd, Double A[], UInt ldA)
{
  F77Int info;
  wrapdpbtrf(upper,static_cast<F77Int>(n), static_cast<F77Int>(kd), A,static_cast<F77Int>(ldA), info);
  return info;
}

inline Int lapack_dgbtrf(UInt n, UInt m, UInt kl, UInt ku, Double A[], UInt ldA, F77Int ipiv[])
{
  F77Int info;
  wrapdgbtrf(static_cast<F77Int>(n), static_cast<F77Int>(m), static_cast<F77Int>(kl), static_cast<F77Int>(ku), A,static_cast<F77Int>(ldA), ipiv,info);
  return info;
}

inline Int lapack_dtbtrs(Bool upper, Bool trans, Bool unitDiag, UInt n, UInt kd, UInt nrhs, const Double A[], UInt ldA, Double B[], UInt ldB)
{
  F77Int info;
  wrapdtbtrs(upper,trans,unitDiag,static_cast<F77Int>(n), static_cast<F77Int>(kd), static_cast<F77Int>(nrhs), A,static_cast<F77Int>(ldA), B,static_cast<F77Int>(ldB), info);
  return info;
}

inline Int lapack_dgbtrs(Bool trans, UInt n, UInt kl, UInt ku, UInt nrhs, const Double A[], UInt ldA, F77Int ipiv[], Double B[], UInt ldB)
{
  F77Int info;
  wrapdgbtrs(trans,static_cast<F77Int>(n), static_cast<F77Int>(kl), static_cast<F77Int>(ku), static_cast<F77Int>(nrhs), A, static_cast<F77Int>(ldA), ipiv, B, static_cast<F77Int>(ldB), info);
  return info;
}

inline Int lapack_dsyev(Bool jobz, Bool upper, UInt n, Double A[], UInt ldA, Double W[])
{
  F77Int    info;
  F77Double tmp;
  wrapdsyev(jobz, upper,static_cast<F77Int>(n), A,static_cast<F77Int>(ldA), W,&tmp,-1,info);
  F77Int lwork = (tmp>0) ? static_cast<F77Int>(tmp) : static_cast<F77Int>(n);
  F77Double *work = new F77Double[lwork];
  wrapdsyev(jobz, upper,static_cast<F77Int>(n), A,static_cast<F77Int>(ldA), W,work,lwork,info);
  delete[] work;
  return info;
}

inline Int lapack_dgeev(Bool jobvl, Bool jobvr, UInt n, Double A[], UInt ldA, Double WR[], Double WI[], Double VL[], UInt ldVL, Double VR[], UInt ldVR)
{
  F77Int    info;
  F77Double tmp;
  wrapdgeev(jobvl,jobvr,static_cast<F77Int>(n), A,static_cast<F77Int>(ldA), WR,WI,VL, static_cast<F77Int>(ldVL), VR, static_cast<F77Int>(ldVR), &tmp,-1,info);
  F77Int lwork = (tmp>0) ? static_cast<F77Int>(tmp) : static_cast<F77Int>(n);
  F77Double *work = new F77Double[lwork];
  wrapdgeev(jobvl,jobvr,static_cast<F77Int>(n), A,static_cast<F77Int>(ldA), WR,WI,VL, static_cast<F77Int>(ldVL), VR, static_cast<F77Int>(ldVR), work,lwork,info);
  delete[] work;
  return info;
}

inline Int lapack_dgesvd(Bool jobu, Bool jobvt, UInt m, UInt n, Double A[], UInt ldA, Double S[], Double U[], UInt ldU, Double VT[], UInt ldvt)
{
  F77Int    info;
  F77Double tmp;
  wrapdgesvd(jobu,jobvt,static_cast<F77Int>(m), static_cast<F77Int>(n), A,static_cast<F77Int>(ldA), S,U, static_cast<F77Int>(ldU), VT, static_cast<F77Int>(ldvt), &tmp,-1,info);
  F77Int lwork = (tmp>0) ? static_cast<F77Int>(tmp) : static_cast<F77Int>(std::max(n,m));
  F77Double *work = new F77Double[lwork];
  wrapdgesvd(jobu,jobvt,static_cast<F77Int>(m), static_cast<F77Int>(n), A,static_cast<F77Int>(ldA), S,U, static_cast<F77Int>(ldU), VT, static_cast<F77Int>(ldvt), work,lwork,info);
  delete[] work;
  return info;
}

inline Int lapack_dgesdd(Bool jobz, UInt m, UInt n, Double A[], UInt ldA, Double S[], Double U[], UInt ldU, Double VT[], UInt ldvt)
{
  F77Int    info;
  F77Double tmp;
  F77Int itmp;
  wrapdgesdd(jobz,static_cast<F77Int>(m), static_cast<F77Int>(n), A,static_cast<F77Int>(ldA), S,U, static_cast<F77Int>(ldU), VT, static_cast<F77Int>(ldvt), &tmp,-1,&itmp,info);
  F77Int lwork = (tmp>0) ? static_cast<F77Int>(tmp) : static_cast<F77Int>(std::max(n,m));
  F77Double *work = new F77Double[lwork];
  F77Int *iwork = new F77Int[8*std::min(m,n)];
  wrapdgesdd(jobz,static_cast<F77Int>(m), static_cast<F77Int>(n), A,static_cast<F77Int>(ldA), S,U, static_cast<F77Int>(ldU), VT, static_cast<F77Int>(ldvt), work,lwork,iwork,info);
  delete[] work;
  delete[] iwork;
  return info;
}

/***********************************************/

#endif /* __GROOPS_LAPACK__ */
