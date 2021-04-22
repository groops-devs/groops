/***********************************************/
/**
* @file kernelTruncation.cpp
*
* @brief Kernel truncated at spherical harmonic degree,
* @see Kernel
*
* @author Andreas Kvas
* @date 2019-03-01
*
*/
/***********************************************/

#include "base/import.h"
#include "classes/kernel/kernel.h"
#include "classes/kernel/kernelTruncation.h"

/***********************************************/

KernelTruncation::KernelTruncation(Config &config)
{
  try
  {
    readConfig(config, "kernel",    _kernel,    Config::MUSTSET, "",  "");
    readConfig(config, "minDegree", _minDegree, Config::DEFAULT, "0", "truncate before minDegree");
    readConfig(config, "maxDegree", _maxDegree, Config::MUSTSET, "", "truncate after maxDegree");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector KernelTruncation::coefficients(const Vector3d &p, UInt degree) const
{
  Vector kn = _kernel->coefficients(p, degree);
  if(_minDegree > 0)
    kn.row(0, _minDegree).setNull();
  if(_maxDegree+2 < kn.rows())
    kn.row(_maxDegree+1, kn.rows()-_maxDegree-1).setNull();
  return kn;
}

/***********************************************/

Vector KernelTruncation::inverseCoefficients(const Vector3d &p, UInt degree, Bool interior) const
{
  Vector kn = _kernel->inverseCoefficients(p, degree, interior);
  if(_minDegree > 0)
    kn.row(0, _minDegree).setNull();
  if(_maxDegree+2 < kn.rows())
    kn.row(_maxDegree+1, kn.rows()-_maxDegree-1).setNull();
  return kn;
}

/***********************************************/
