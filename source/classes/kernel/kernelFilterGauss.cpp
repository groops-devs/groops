/***********************************************/
/**
* @file kernelFilterGauss.cpp
*
* @brief Kernel smoothed by an Gaussian filter.
* @see Kernel
*
* @author Torsten Mayer-Guerr
* @date 2007-10-10
*
*/
/***********************************************/

#include "base/import.h"
#include "classes/kernel/kernel.h"
#include "classes/kernel/kernelFilterGauss.h"

/***********************************************/

KernelFilterGauss::KernelFilterGauss(Config &config)
{
  try
  {
    readConfig(config, "kernel", kernel, Config::MUSTSET, "",    "");
    readConfig(config, "radius", radius, Config::MUSTSET, "500", "filter radius [km]");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector KernelFilterGauss::computeFilterCoefficients(UInt degree) const
{
  try
  {
    Double b = log(2.)/(1-cos(1000*radius/DEFAULT_R));

    Vector wn(degree+1);
    wn(0) = 1;
    wn(1) = (1+exp(-2*b))/(1-exp(-2*b)) - 1./b;
    for(UInt n=2; n<=degree; n++)
    {
      wn(n) = -(2*n-1.)/b * wn(n-1) + wn(n-2);
      if(wn(n)<1e-7)
        break;
    }

    return wn;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector KernelFilterGauss::coefficients(const Vector3d &p, UInt degree) const
{
  Vector kn = kernel->coefficients(p,degree);
  if(filterCoeff.size()<kn.size())
    filterCoeff = computeFilterCoefficients(kn.size()-1);
  for(UInt n=0; n<kn.size(); n++)
    kn(n) *= filterCoeff(n);
  return kn;
}

/***********************************************/

Vector KernelFilterGauss::inverseCoefficients(const Vector3d &p, UInt degree, Bool interior) const
{
  Vector kn = kernel->inverseCoefficients(p, degree, interior);
  if(filterCoeff.size()<kn.size())
    filterCoeff = computeFilterCoefficients(kn.size()-1);
  for(UInt n=0; n<kn.size(); n++)
    kn(n) *= filterCoeff(n);
  return kn;
}

/***********************************************/
