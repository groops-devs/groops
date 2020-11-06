/***********************************************/
/**
* @file kernelBlackmanLowPass.cpp
*
* @brief Kernel smoothed by Blackman low-pass filter
* @see Kernel
*
* @author Andreas Kvas
* @date 2019-03-01
*
*/
/***********************************************/

#include "base/import.h"
#include "classes/kernel/kernel.h"
#include "classes/kernel/kernelBlackmanLowPass.h"

/***********************************************/

KernelBlackmanLowPass::KernelBlackmanLowPass(Config &config)
{
  try
  {
    readConfig(config, "kernel",                kernel, Config::MUSTSET, "", "");
    readConfig(config, "startDegreeTransition", n1,     Config::MUSTSET, "", "minimum degree in transition band");
    readConfig(config, "stopDegreeTransition",  n2,     Config::MUSTSET, "", "maximum degree in transition band");
    if(isCreateSchema(config)) return;

    if(n1>n2)
      throw(Exception("Transition band limits must be given in increasing order."));

  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector KernelBlackmanLowPass::computeFilterCoefficients(UInt degree) const
{
  try
  {
    Vector wn(degree+1);
    for(UInt n = 0; n < std::min(n1, degree+1); n++)  // pass band
      wn(n) = 1.0;

    for(UInt n = n1; n <= std::min(n2, degree); n++)
    {
      wn(n) = std::pow(0.42 + 0.5 * std::cos(PI*static_cast<Double>(n-n1)/static_cast<Double>(n2-n1)) +
      0.08 * std::cos(2*PI*static_cast<Double>(n-n1)/static_cast<Double>(n2-n1)), 1);
    }

    return wn;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector KernelBlackmanLowPass::coefficients(const Vector3d &p, UInt degree) const
{
  Vector kn = kernel->coefficients(p,degree);
  if(filterCoeff.size()<kn.size())
    filterCoeff = computeFilterCoefficients(kn.size()-1);
  for(UInt n=0; n<kn.size(); n++)
    kn(n) *= filterCoeff(n);
  return kn;
}

/***********************************************/

Vector KernelBlackmanLowPass::inverseCoefficients(const Vector3d &p, UInt degree, Bool interior) const
{
  Vector kn = kernel->inverseCoefficients(p, degree, interior);
  if(filterCoeff.size()<kn.size())
    filterCoeff = computeFilterCoefficients(kn.size()-1);
  for(UInt n=0; n<kn.size(); n++)
    kn(n) *= filterCoeff(n);
  return kn;
}

/***********************************************/
