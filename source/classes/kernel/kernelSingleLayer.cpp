/***********************************************/
/**
* @file kernelSingleLayer.cpp
*
* @brief Single layer density.
* @see Kernel
*
* @author Torsten Mayer-Guerr
* @date 2003-09-20
*
*/
/***********************************************/

#include "base/import.h"
#include "base/legendrePolynomial.h"
#include "files/fileMatrix.h"
#include "classes/gravityfield/gravityfield.h"
#include "classes/kernel/kernel.h"
#include "classes/kernel/kernelSingleLayer.h"

/***********************************************/

KernelSingleLayer::KernelSingleLayer(Config &config)
{
  try
  {
    FileName loveNumberName;

    readConfig(config, "inputfileLoadingLoveNumber", loveNumberName, Config::OPTIONAL,  "{groopsDataDir}/loading/loadLoveNumbers_Gegout97.txt", "");
    if(isCreateSchema(config)) return;

    if(!loveNumberName.empty())
      readFileMatrix(loveNumberName, kn);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector KernelSingleLayer::coefficients(Vector3d const &q, UInt degree) const
{
  try
  {
    if((degree==INFINITYDEGREE) && (kn.rows()>0))
      degree = kn.rows()-1;
    if(degree==INFINITYDEGREE)
      throw(Exception("INFINITYDEGREE requested"));

    Vector  coeff(degree+1);
    Double *cp = coeff.field();
    Double  factor = 4*PI*GRAVITATIONALCONSTANT*q.r();
    for(UInt n=0; n<std::min(degree+1, kn.size()); n++)
      *cp++ = factor * (1.+kn(n)) / (2.*n+1.);
    for(UInt n=std::min(degree+1, kn.size()); n<=degree; n++)
      *cp++ = factor / (2.*n+1.);

    return coeff;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector KernelSingleLayer::inverseCoefficients(Vector3d const &p, UInt degree, Bool interior) const
{
  try
  {
    if((degree==INFINITYDEGREE) && (kn.rows()>0))
      degree = kn.rows()-1;
    if(degree==INFINITYDEGREE)
      throw(Exception("INFINITYDEGREE requested"));

    if(interior)
      throw(Exception("interior not implemented"));

    Vector  coeff(degree+1);
    Double *cp = coeff.field();
    Double  factor = 1./(4.*PI*GRAVITATIONALCONSTANT*p.r());
    for(UInt n=0; n<std::min(degree+1, kn.size()); n++)
      *cp++ = factor * (2.*n+1.) / (1.+kn(n));
    for(UInt n=std::min(degree+1, kn.size()); n<=degree; n++)
      *cp++ = factor * (2.*n+1.);

    return coeff;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double KernelSingleLayer::kernel(Vector3d const &p, Vector3d const &q) const
{
  try
  {
    Double K = 4*PI*GRAVITATIONALCONSTANT/(p-q).r();

    // add love load numbers
    if(kn.size())
    {
      Vector  coeff(kn.rows());
      Double *cp = coeff.field();
      const Double factor = 4*PI*GRAVITATIONALCONSTANT/q.r();
      for(UInt n=0; n<kn.rows(); n++)
        *cp++ = factor*kn(n)/(2.*n+1.);
      K += Kernel::kernel(p, q, coeff);
    }

    return K;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double KernelSingleLayer::radialDerivative(Vector3d const &p, Vector3d const &q) const
{
  try
  {
    const Double r = p.r();
    const Double R = q.r();
    const Double l = (p-q).r();
    const Double cos_psi = (R*R+r*r-l*l)/(2*R*r);
    Double dKdr = -4*PI*GRAVITATIONALCONSTANT/(l*l*l) * (r-R*cos_psi);

    if(kn.rows())
    {
      Vector  coeff(kn.rows());
      Double *cp = coeff.field();
      const Double factor = 4*PI*GRAVITATIONALCONSTANT/q.r();
      for(UInt n=0; n<kn.rows(); n++)
        *cp++ = factor*kn(n)/(2.*n+1.);
      dKdr += Kernel::radialDerivative(p, q, coeff);
    }

    return dKdr;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector3d KernelSingleLayer::gradient(Vector3d const &p, Vector3d const &q) const
{
  try
  {
    const Double l = (p-q).r();
    Vector3d g = -4*PI*GRAVITATIONALCONSTANT/(l*l*l) * (p-q);

    if(kn.rows())
    {
      Vector  coeff(kn.rows());
      Double *cp = coeff.field();
      const Double factor = 4*PI*GRAVITATIONALCONSTANT/q.r();
      for(UInt n=0; n<kn.rows(); n++)
        *cp++ = factor*kn(n)/(2.*n+1.);
      g += Kernel::gradient(p, q, coeff);
    }

    return g;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double KernelSingleLayer::inverseKernel(Vector3d const &p, Vector3d const &q, const Kernel &kernel) const
{
  try
  {
    Double f = -1./(4.*PI*GRAVITATIONALCONSTANT) * (2*kernel.radialDerivative(p,q) + kernel.kernel(p,q)/p.r());

    if(kn.rows())
    {
      Double r  = p.r();
      Double R  = q.r();
      Double t  = inner(p, q)/r/R; // t = cos(psi)
      Vector k2 = kernel.coefficients(p, kn.rows()-1);
      Vector k1(kn.rows());
      Double f1 = R/r;
      Double f2 = R/r;
      Double factor = 1./(4.*PI*GRAVITATIONALCONSTANT*p.r());
      Double       *p1 = k1.field();
      const Double *p2 = k2.field();
      for(UInt n=0; n<kn.rows(); n++)
      {
        // k1(n) *= (R/r)^(n+1) * sqrt(2n+1) * k2(n));
        *p1++ = f1 * factor * (2.*n+1.) * (1/(1.+kn(n))-1) * sqrt(2.*n+1.) * *p2++;
        f1  *= f2;
      }
      f += LegendrePolynomial::sum(t, k1, kn.rows()-1);
    }

    return f;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double KernelSingleLayer::inverseKernel(const Time &time, const Vector3d &p, const GravityfieldBase &field) const
{
  try
  {
    Double f = -1./(4.*PI*GRAVITATIONALCONSTANT) * (2*field.radialGradient(time,p) + field.potential(time,p)/p.r());

    if(kn.rows())
    {
      Vector  coeff(kn.rows());
      Double *cp = coeff.field();
      Double  factor = 1./(4.*PI*GRAVITATIONALCONSTANT*p.r());
      for(UInt n=0; n<kn.rows(); n++)
        *cp++ = factor * (2.*n+1.) * (1/(1.+kn(n))-1);
      // Convolution with the kernel
      SphericalHarmonics harmonics = field.sphericalHarmonics(time, coeff.rows()-1);
      f += inner(coeff, harmonics.Yn(p));
    }

    return f;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
