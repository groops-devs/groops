/***********************************************/
/**
* @file kernelDeformation.cpp
*
* @brief Radial deformation by loading.
* @see Kernel
*
* @author Torsten Mayer-Guerr
* @date 2009-07-29
*
*/
/***********************************************/

#include "base/import.h"
#include "base/planets.h"
#include "files/fileMatrix.h"
#include "classes/kernel/kernel.h"
#include "classes/kernel/kernelDeformation.h"

/***********************************************/

KernelDeformation::KernelDeformation(Config &config)
{
  try
  {
    FileName deformationName, potentialName;

    readConfig(config, "inputfileDeformationLoadLoveNumber", deformationName, Config::MUSTSET,   "{groopsDataDir}/loading/deformationLoveNumbers_CM_Gegout97.txt", "");
    readConfig(config, "inputfilePotentialLoadLoveNumber",   potentialName,   Config::OPTIONAL,  "{groopsDataDir}/loading/loadLoveNumbers_Gegout97.txt", "if full potential is given and not only loading potential");
    if(isCreateSchema(config)) return;

    // deformation load love numbers (hn,ln)
    // -------------------------------------
    readFileMatrix(deformationName, love);

    // Compute loading potenial, if full potential is given
    // ----------------------------------------------------
    if(!potentialName.empty())
    {
      Vector kn;
      readFileMatrix(potentialName, kn);

      for(UInt n=2; n<std::min(kn.rows(), love.rows()); n++)
        love(n,0) /= (1.+kn(n));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector KernelDeformation::coefficients(Vector3d const &q, UInt degree) const
{
  if(degree==INFINITYDEGREE)
    degree = love.rows()-1;

  Vector  k(degree+1);
  Double *kn = k.field();
  Double  factor = Planets::normalGravity(q);
  for(UInt n=0; n<std::min(degree+1, love.size()); n++)
    *kn++ = ((love(n,0) == 0.) ? (0.0) : (factor / love(n,0)));

  return k;
}

/***********************************************/

Vector KernelDeformation::inverseCoefficients(Vector3d const &p, UInt degree, Bool /*interior*/) const
{
  if(degree==INFINITYDEGREE)
    degree = love.rows()-1;

  Vector  k(degree+1);
  Double *kn = k.field();
  Double  factor = 1./Planets::normalGravity(p);
  for(UInt n=0; n<std::min(degree+1, love.size()); n++)
    *kn++ = factor * love(n,0);

  return k;
}

/***********************************************/
/***********************************************/
