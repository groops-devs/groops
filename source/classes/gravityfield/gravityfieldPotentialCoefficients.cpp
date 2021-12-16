/***********************************************/
/**
* @file gravityfieldPotentialCoefficients.cpp
*
* @brief Potential coefficients (SphericalHarmonics).
* All classes fo this type are added together before evaluation to speed up the computation.
* @see Gravityfield
*
* @author Torsten Mayer-Guerr
* @date 2001-08-25
*
*/
/***********************************************/

#include "base/import.h"
#include "base/sphericalHarmonics.h"
#include "config/config.h"
#include "files/fileSphericalHarmonics.h"
#include "classes/kernel/kernel.h"
#include "classes/gravityfield/gravityfield.h"
#include "classes/gravityfield/gravityfieldPotentialCoefficients.h"

/***********************************************/

GravityfieldPotentialCoefficients::GravityfieldPotentialCoefficients(Config &config)
{
  try
  {
    FileName fileName;
    Double   factor;
    UInt     minDegree, maxDegree = INFINITYDEGREE;
    Bool     zeroVariance;

    readConfig(config, "inputfilePotentialCoefficients", fileName, Config::MUSTSET, "{groopsDataDir}/potential/", "");
    readConfig(config, "minDegree",       minDegree,    Config::DEFAULT,  "0",   "");
    readConfig(config, "maxDegree",       maxDegree,    Config::OPTIONAL, "",    "");
    readConfig(config, "factor",          factor,       Config::DEFAULT,  "1.0", "the result is multiplied by this factor, set -1 to subtract the field");
    readConfig(config, "setSigmasToZero", zeroVariance, Config::DEFAULT,  "0",   "set variances to zero, should be used by adding back reference fields");
    if(isCreateSchema(config)) return;

    readFileSphericalHarmonics(fileName, harmonics);
    if(zeroVariance)
      harmonics.sigma2cnm() = harmonics.sigma2snm() = Matrix();
    harmonics = factor * harmonics.get(maxDegree, minDegree);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double GravityfieldPotentialCoefficients::potential(const Time &/*time*/, const Vector3d &point) const
{
  return harmonics.potential(point);
}

/***********************************************/

Double GravityfieldPotentialCoefficients::radialGradient(const Time &/*time*/, const Vector3d &point) const
{
  return harmonics.radialGradient(point);
}

/***********************************************/

Double GravityfieldPotentialCoefficients::field(const Time &/*time*/, const Vector3d &point, const Kernel &kernel) const
{
  // Convolution with kernel
  return inner(kernel.inverseCoefficients(point, harmonics.maxDegree(), harmonics.isInterior()), harmonics.Yn(point, harmonics.maxDegree()));
}

/***********************************************/

Vector3d GravityfieldPotentialCoefficients::gravity(const Time &/*time*/, const Vector3d &point) const
{
  return harmonics.gravity(point);
}

/***********************************************/

Tensor3d GravityfieldPotentialCoefficients::gravityGradient(const Time &/*time*/, const Vector3d &point) const
{
  return harmonics.gravityGradient(point);
}

/***********************************************/

Vector3d GravityfieldPotentialCoefficients::deformation(const Time &/*time*/, const Vector3d &point, Double gravity, const Vector &hn, const Vector &ln) const
{
  return harmonics.deformation(point, gravity, hn, ln);
}

/***********************************************/

void GravityfieldPotentialCoefficients::deformation(const std::vector<Time> &time, const std::vector<Vector3d> &point, const std::vector<Double> &gravity,
                                                         const Vector &hn, const Vector &ln, std::vector< std::vector<Vector3d> > &disp) const
{
  for(UInt k=0; k<point.size(); k++)
  {
    Vector3d d = harmonics.deformation(point.at(k), gravity.at(k), hn, ln);
    for(UInt i=0; i<time.size(); i++)
      disp.at(k).at(i) += d;
  }
}

/***********************************************/

void GravityfieldPotentialCoefficients::variance(const Time &/*time*/, const std::vector<Vector3d> &point, const Kernel &kernel, Matrix &D) const
{
  try
  {
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

Double GravityfieldPotentialCoefficients::variance(const Time &/*time*/, const Vector3d &point, const Kernel &kernel) const
{
  try
  {
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
        sigma2  += pow(coeff(n)*Cnm(n,m),2) * sigma2cnm(n,m)
                +  pow(coeff(n)*Snm(n,m),2) * sigma2snm(n,m);
    return sigma2;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

SphericalHarmonics GravityfieldPotentialCoefficients::sphericalHarmonics(const Time &/*time*/, UInt maxDegree, UInt minDegree, Double GM, Double R) const
{
  try
  {
    return harmonics.get(maxDegree, minDegree, GM, R);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
