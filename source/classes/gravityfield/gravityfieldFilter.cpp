/***********************************************/
/**
* @file gravityfieldFilter.cpp
*
* @brief Filtered spherical harmonics.
* @see Gravityfield
*
* @author Torsten Mayer-Guerr
* @date 2020-07-20
*
*/
/***********************************************/

#include "base/import.h"
#include "base/sphericalHarmonics.h"
#include "config/config.h"
#include "classes/kernel/kernel.h"
#include "classes/gravityfield/gravityfield.h"
#include "classes/gravityfield/gravityfieldFilter.h"

/***********************************************/

GravityfieldFilter::GravityfieldFilter(Config &config)
{
  try
  {
    readConfig(config, "gravityfield", gravityfield, Config::MUSTSET, "", "");
    readConfig(config, "filter",       filter,       Config::MUSTSET, "", "");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double GravityfieldFilter::potential(const Time &time, const Vector3d &point) const
{
  return sphericalHarmonics(time).potential(point);
}

/***********************************************/

Double GravityfieldFilter::radialGradient(const Time &time, const Vector3d &point) const
{
  return sphericalHarmonics(time).radialGradient(point);
}

/***********************************************/

Double GravityfieldFilter::field(const Time &time, const Vector3d &point, const Kernel &kernel) const
{
  SphericalHarmonics harmonics = sphericalHarmonics(time);
  return inner(kernel.inverseCoefficients(point, harmonics.maxDegree(), harmonics.isInterior()), harmonics.Yn(point, harmonics.maxDegree()));
}

/***********************************************/

Vector3d GravityfieldFilter::gravity(const Time &time, const Vector3d &point) const
{
  return sphericalHarmonics(time).gravity(point);
}

/***********************************************/

Tensor3d GravityfieldFilter::gravityGradient(const Time &time, const Vector3d &point) const
{
  return sphericalHarmonics(time).gravityGradient(point);
}

/***********************************************/

Vector3d GravityfieldFilter::deformation(const Time &time, const Vector3d &point, Double gravity, const Vector &hn, const Vector &ln) const
{
  return sphericalHarmonics(time).deformation(point, gravity, hn, ln);
}

/***********************************************/

void GravityfieldFilter::deformation(const std::vector<Time> &time, const std::vector<Vector3d> &point, const std::vector<Double> &gravity,
                                     const Vector &hn, const Vector &ln, std::vector<std::vector<Vector3d>> &disp) const
{
  for(UInt i=0; i<time.size(); i++)
  {
    SphericalHarmonics harmonics = sphericalHarmonics(time.at(i));
    for(UInt k=0; k<point.size(); k++)
      disp.at(k).at(i) += harmonics.deformation(point.at(k), gravity.at(k), hn, ln);
  }
}

/***********************************************/

SphericalHarmonics GravityfieldFilter::sphericalHarmonics(const Time &time, UInt maxDegree, UInt minDegree, Double GM, Double R) const
{
  return filter->filter(gravityfield->sphericalHarmonics(time, maxDegree, minDegree, GM, R));
}

/***********************************************/

void GravityfieldFilter::variance(const Time &time, const std::vector<Vector3d> &point, const Kernel &kernel, Matrix &D) const
{
  try
  {
    SphericalHarmonics harmonics = sphericalHarmonics(time);
    if(!harmonics.sigma2cnm().size())
      return;

    UInt   maxDegree = harmonics.maxDegree();
    Double GM = harmonics.GM();
    Double R  = harmonics.R();

    Matrix Cnm_i, Snm_i, Cnm_k, Snm_k;
    Matrix sigma2cnm = harmonics.sigma2cnm();
    Matrix sigma2snm = harmonics.sigma2snm();

    for(UInt i=0; i<point.size(); i++)
    {
      Vector coeff_i = GM/R * kernel.inverseCoefficients(point.at(i), maxDegree, harmonics.isInterior());
      SphericalHarmonics::CnmSnm(1/R * point.at(i), maxDegree, Cnm_i, Snm_i, harmonics.isInterior());

      for(UInt k=i; k<point.size(); k++)
      {
        Vector coeff_k = GM/R * kernel.inverseCoefficients(point.at(k), maxDegree, harmonics.isInterior());
        SphericalHarmonics::CnmSnm(1/R * point.at(k), maxDegree, Cnm_k, Snm_k, harmonics.isInterior());

        for(UInt n=0; n<=maxDegree; n++)
          for(UInt m=0; m<=n; m++)
            D(i,k) += coeff_i(n) * Cnm_i(n,m) * sigma2cnm(n,m) * Cnm_k(n,m) * coeff_k(n)
                   +  coeff_i(n) * Snm_i(n,m) * sigma2snm(n,m) * Snm_k(n,m) * coeff_k(n);
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double GravityfieldFilter::variance(const Time &time, const Vector3d &point, const Kernel &kernel) const
{
  try
  {
    SphericalHarmonics harmonics = sphericalHarmonics(time);
    if(!harmonics.sigma2cnm().size())
      return 0.0;

    UInt   maxDegree = harmonics.maxDegree();
    Double GM = harmonics.GM();
    Double R  = harmonics.R();
    Vector coeff = GM/R * kernel.inverseCoefficients(point, maxDegree, harmonics.isInterior());
    Matrix Cnm, Snm;
    SphericalHarmonics::CnmSnm(1/R * point, maxDegree, Cnm, Snm, harmonics.isInterior());

    Double sigma2 = 0;
    const Matrix sigma2cnm = harmonics.sigma2cnm();
    const Matrix sigma2snm = harmonics.sigma2snm();
    for(UInt n=0; n<=maxDegree; n++)
      for(UInt m=0; m<=n; m++)
        sigma2  += std::pow(coeff(n)*Cnm(n,m),2) * sigma2cnm(n,m)
                +  std::pow(coeff(n)*Snm(n,m),2) * sigma2snm(n,m);
    return sigma2;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
