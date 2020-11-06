/***********************************************/
/**
* @file sphericalHarmonicsFilterGauss.h
*
* @brief Isotropic gaussian filter.
* @see SphericalHarmonicsFilter
*
* @author Torsten Mayer-Guerr
* @date 2008-08-06
*
*/
/***********************************************/

#ifndef __GROOPS_SPHERICALHARMONICSFILTERGAUSS__
#define __GROOPS_SPHERICALHARMONICSFILTERGAUSS__

// Latex documentation
#ifdef DOCSTRING_SphericalHarmonicsFilter
static const char *docstringSphericalHarmonicsFilterGauss = R"(
\subsection{Gauss}
Filtering the spherical harmonics expansion with a Gaussian filter.
\config{radius} gives the filter radius on the Earth surface in km.
)";
#endif

/***********************************************/

#include "classes/sphericalHarmonicsFilter/sphericalHarmonicsFilter.h"

/***** CLASS ***********************************/

/** @brief Isotropic gaussian filter.
* @ingroup sphericalHarmonicsFilterGroup
* @see SphericalHarmonicsFilter */
class SphericalHarmonicsFilterGauss : public SphericalHarmonicsFilterBase
{
  mutable Vector filterCoeff;
  Double  radius;

  Vector computeFilterCoefficients(UInt degree) const;

public:
  SphericalHarmonicsFilterGauss(Config &config);

  SphericalHarmonics filter(const SphericalHarmonics &harm) const;
};

/***********************************************/

inline SphericalHarmonicsFilterGauss::SphericalHarmonicsFilterGauss(Config &config)
{
  try
  {
    readConfig(config, "radius", radius, Config::MUSTSET, "500", "filter radius [km]");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Vector SphericalHarmonicsFilterGauss::computeFilterCoefficients(UInt degree) const
{
  try
  {
    const Double b = std::log(2.)/(1-std::cos(1000*radius/DEFAULT_R));

    Vector wn(degree+1);
    wn(0) = 1;
    wn(1) = (1+std::exp(-2*b))/(1-std::exp(-2*b)) - 1./b;
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

inline SphericalHarmonics SphericalHarmonicsFilterGauss::filter(const SphericalHarmonics &harm) const
{
  try
  {
    if(harm.maxDegree()+1 > filterCoeff.size())
      filterCoeff = computeFilterCoefficients(harm.maxDegree());

    Matrix cnm       = harm.cnm();
    Matrix snm       = harm.snm();
    Matrix sigma2cnm = harm.sigma2cnm();
    Matrix sigma2snm = harm.sigma2snm();
    for(UInt n=0; n<cnm.rows(); n++)
    {
      cnm.slice(n,0,1,n+1)       *= filterCoeff(n);
      snm.slice(n,0,1,n+1)       *= filterCoeff(n);
      if(sigma2cnm.size()) sigma2cnm.slice(n,0,1,n+1) *= pow(filterCoeff(n),2);
      if(sigma2snm.size()) sigma2snm.slice(n,0,1,n+1) *= pow(filterCoeff(n),2);
    }
    return SphericalHarmonics(harm.GM(), harm.R(), cnm, snm, sigma2cnm, sigma2snm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
