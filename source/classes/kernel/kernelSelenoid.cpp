/***********************************************/
/**
* @file kernelSelenoid.cpp
*
* @brief Integral kernel of selenoid computation.
* (= Poisson Kern * gamma).
* @see Kernel
*
* @author Beate Klinger
* @date 2013-xx-xx
*
*/
/***********************************************/

#include "base/import.h"
#include "classes/gravityfield/gravityfield.h"
#include "classes/kernel/kernel.h"
#include "classes/kernel/kernelSelenoid.h"

/***********************************************/

Vector KernelSelenoid::coefficients(Vector3d const &/*p*/, UInt degree) const
{
  try
  {
    if(degree==INFINITYDEGREE)
      throw(Exception("INFINITYDEGREE requested"));

    Vector  k(degree+1);
    Double *kn = k.field();
    const Double gamma = DEFAULT_GM/(pow(DEFAULT_R,2));
    for(UInt n=0; n<=degree; n++)
      *kn++ = gamma;

    return k;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector KernelSelenoid::inverseCoefficients(Vector3d const &/*p*/, UInt degree, Bool /*interior*/) const
{
  try
  {
    if(degree==INFINITYDEGREE)
      throw(Exception("INFINITYDEGREE requested"));

    Vector  k(degree+1);
    Double *kn = k.field();
    const Double gamma = DEFAULT_GM/(pow(DEFAULT_R,2));
    for(UInt n=0; n<=degree; n++)
      *kn++ = 1./gamma;

    return k;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double KernelSelenoid::kernel(Vector3d const &p, Vector3d const &q) const
{
  Vector3d diff  = p-q;
  Double r       = p.r();
  Double R       = q.r();
  Double l       = diff.r();

  const Double gamma = DEFAULT_GM/(pow(DEFAULT_R,2));

  return gamma * R*(r*r-R*R)/(l*l*l);
}

/***********************************************/

Double KernelSelenoid::radialDerivative(Vector3d const &p, Vector3d const &q) const
{
  Vector3d diff  = p-q;
  Double r       = p.r();
  Double R       = q.r();
  Double l       = diff.r();
  Double cos_psi = (R*R+r*r-l*l)/(2*R*r);
  Double dl_dr   = (r-R*cos_psi)/l;

  const Double gamma = DEFAULT_GM/(pow(DEFAULT_R,2));

  return gamma * R*(2*r*l-(r*r-R*R)*3*dl_dr)/(l*l*l*l);
}

/***********************************************/

Vector3d KernelSelenoid::gradient(Vector3d const &p, Vector3d const &q) const
{
  Vector3d diff  = p-q;
  Double r       = p.r();
  Double l       = diff.r();
  Double R       = q.r();

  const Double gamma = DEFAULT_GM/(pow(DEFAULT_R,2));

  return gamma * (2/(l*l*l)*p - 3*(r*r-R*R)/(l*l*l*l*l)*diff);
}

/***********************************************/

Tensor3d KernelSelenoid::gradientGradient(const Vector3d &p, const Vector3d &q) const
{
  Tensor3d tns;

  // Vorausberechnungen
  Vector3d diff = p-q;
  Double r      = p.r();
  Double R      = q.r();
  Double l      = diff.r();
  Double l2     = l*l;
  Double l3     = l2*l;
  Double l5     = l3*l2;
  Double l7     = l5*l2;
  Double term   = r*r-R*R;

  const Double gamma = DEFAULT_GM/(pow(DEFAULT_R,2));

  tns.xx() =  2.0/l3 - 12.0*p.x()*diff.x()/l5 + 15*term*diff.x()*diff.x()/l7 - 3.0*term/l5;
  tns.yy() =  2.0/l3 - 12.0*p.y()*diff.y()/l5 + 15*term*diff.y()*diff.y()/l7 - 3.0*term/l5;
  tns.zz() =  2.0/l3 - 12.0*p.z()*diff.z()/l5 + 15*term*diff.z()*diff.z()/l7 - 3.0*term/l5;
  tns.xy() = -6.0*(p.x()*diff.y()+p.y()*diff.x())/l5 + 15.0*term*diff.x()*diff.y()/l7;
  tns.xz() = -6.0*(p.x()*diff.z()+p.z()*diff.x())/l5 + 15.0*term*diff.x()*diff.z()/l7;
  tns.yz() = -6.0*(p.y()*diff.z()+p.z()*diff.y())/l5 + 15.0*term*diff.y()*diff.z()/l7;

  return gamma * tns;
}

/***********************************************/

Double KernelSelenoid::inverseKernel(Vector3d const &p, Vector3d const &q, const Kernel &kernel) const
{
  // potential = K
  const Double gamma = DEFAULT_GM/(pow(DEFAULT_R,2));

  return (1./gamma) * kernel.kernel(p, q);
}

/***********************************************/

Double KernelSelenoid::inverseKernel(const Time &time, const Vector3d &p, const GravityfieldBase &field) const
{
  const Double gamma = DEFAULT_GM/(pow(DEFAULT_R,2));

  return (1./gamma) * field.potential(time, p);
}

/***********************************************/
