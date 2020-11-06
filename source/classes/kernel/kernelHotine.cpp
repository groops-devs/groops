/***********************************************/
/**
* @file kernelHotine.cpp
*
* @brief Hotine kernel (gravity disturbances).
* @see Kernel
*
* @author Torsten Mayer-Guerr
* @date 2003-09-20
*
*/
/***********************************************/


#include "base/import.h"
#include "classes/gravityfield/gravityfield.h"
#include "classes/kernel/kernel.h"
#include "classes/kernel/kernelHotine.h"

/***********************************************/

Vector KernelHotine::coefficients(Vector3d const &q, UInt degree) const
{
  if(degree==INFINITYDEGREE)
    throw(Exception("In KernelHotine::coefficients: INFINITYDEGREE requested"));

  Double  R = q.r();
  Vector  k(degree+1);
  Double *kn = k.field();

  for(UInt n=0; n<=degree; n++)
    *kn++ = R/(n+1);

  return k;
}

/***********************************************/

Vector KernelHotine::inverseCoefficients(Vector3d const &p, UInt degree, Bool interior) const
{
  if(degree==INFINITYDEGREE)
    throw(Exception("In KernelHotine::inverseCoefficients: INFINITYDEGREE requested"));

  Double  r = p.r();
  Vector  k(degree+1);
  Double *kn = k.field();

  if(interior)
    for(UInt n=0; n<=degree; n++)
      *kn++ = -(n/r);
  else
    for(UInt n=0; n<=degree; n++)
      *kn++ = (n+1)/r;

  return k;
}

/***********************************************/

Double KernelHotine::kernel(Vector3d const &p, Vector3d const &q) const
{
  Vector3d diff  = p-q;
  Double r       = p.r();
  Double R       = q.r();
  Double l       = diff.r();
  Double cos_psi = (R*R+r*r-l*l)/(2*R*r);

  return 2*R*R/l-R*log((l+R-r*cos_psi)/(r-r*cos_psi));
}

/***********************************************/

Double KernelHotine::radialDerivative(Vector3d const &p, Vector3d const &q) const
{
  Vector3d diff  = p-q;
  Double r       = p.r();
  Double R       = q.r();
  Double l       = diff.r();

  return -R*R*(r*r-R*R)/(l*l*l*r);
}

/***********************************************/

Double KernelHotine::inverseKernel(Vector3d const &p, Vector3d const &q, const Kernel &kernel) const
{
  // gravity disturbance = -dK/dr
  return -kernel.radialDerivative(p,q);
}

/***********************************************/

Double KernelHotine::inverseKernel(const Time &time, const Vector3d &p, const GravityfieldBase &field) const
{
  // gravity disturbance = -dK/dr
  return -field.radialGradient(time, p);
}

/***********************************************/
