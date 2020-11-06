/***********************************************/
/**
* @file kernelCoefficients.cpp
*
* @brief Kernel from coefficients.
* Given as series of Legendre polynomials.
* @see Kernel
*
* @author Torsten Mayer-Guerr
* @date 2003-09-20
*
*/
/***********************************************/

#include "base/import.h"
#include "files/fileMatrix.h"
#include "classes/kernel/kernel.h"
#include "classes/kernel/kernelCoefficients.h"

/***********************************************/

KernelCoefficients::KernelCoefficients(Config &config)
{
  try
  {
    FileName knName;

    readConfig(config, "inputfileCoefficients", knName, Config::MUSTSET, "", "");
    if(isCreateSchema(config)) return;

    readFileMatrix(knName, kn);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector KernelCoefficients::coefficients(const Vector3d &/*p*/, UInt degree) const
{
  if(degree==INFINITYDEGREE)
    degree = kn.rows()-1;

  Vector  k(degree+1);
  Double *kp = k.field();
  const Double *kf = kn.field();
  for(UInt n=0; n<=std::min(degree, kn.rows()-1); n++)
    *kp++ = *kf++;

  return k;
}

/***********************************************/

Vector KernelCoefficients::inverseCoefficients(const Vector3d &/*p*/, UInt degree, Bool /*interior*/) const
{
  if(degree==INFINITYDEGREE)
    degree = kn.rows()-1;

  Vector  k(degree+1);
  Double *kp = k.field();
  const Double *kf = kn.field();
  for(UInt n=0; n<=std::min(degree, kn.rows()-1); n++)
  {
    *kp++ = ((*kf == 0.) ? (0.0) : (1.0/ *kf));
    kf++;
  }

  return k;
}

/***********************************************/
