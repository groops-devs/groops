/***********************************************/
/**
* @file kernelRadialGradient.cpp
*
* @brief 2nd radial derivative (gravity gradients).
* @see Kernel
*
* @author Torsten Mayer-Guerr
* @date 2010-04-27
*
*/
/***********************************************/

#include "base/import.h"
#include "classes/kernel/kernel.h"
#include "classes/kernel/kernelRadialGradient.h"

/***********************************************/

Vector KernelRadialGradient::coefficients(Vector3d const &q, UInt degree) const
{
  if(degree==INFINITYDEGREE)
    throw(Exception("In KernelRadialGradient::coefficients: INFINITYDEGREE requested"));

  Double  R2 = pow(q.r(),2);
  Vector  k(degree+1);
  Double *kn = k.field();

  for(UInt n=0; n<=degree; n++)
    *kn++ = R2/((n+1)*(n+2));

  return k;
}

/***********************************************/

Vector KernelRadialGradient::inverseCoefficients(Vector3d const &p, UInt degree, Bool interior) const
{
  if(degree==INFINITYDEGREE)
    throw(Exception("In KernelRadialGradient::inverseCoefficients: INFINITYDEGREE requested"));

  Double  r2 = pow(p.r(),2);
  Vector  k(degree+1);
  Double *kn = k.field();

  if(interior)
    for(UInt n=2; n<=degree; n++)
      *kn++ = n*(n-1.)/r2;
  else
    for(UInt n=0; n<=degree; n++)
      *kn++ = (n+1)*(n+2)/r2;

  return k;
}

/***********************************************/
