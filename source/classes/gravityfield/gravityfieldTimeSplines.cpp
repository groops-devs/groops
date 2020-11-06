/***********************************************/
/**
* @file gravityfieldTimeSplines.cpp
*
* @brief TimeSplines.
* @see Gravityfield
*
* @author Torsten Mayer-Guerr
* @date 2004-04-14
*
*/
/***********************************************/

#include "base/import.h"
#include "base/sphericalHarmonics.h"
#include "config/config.h"
#include "files/fileTimeSplinesGravityfield.h"
#include "classes/kernel/kernel.h"
#include "classes/gravityfield/gravityfield.h"
#include "classes/gravityfield/gravityfieldTimeSplines.h"

/***********************************************/

GravityfieldTimeSplines::GravityfieldTimeSplines(Config &config)
{
  try
  {
    maxDegree = INFINITYDEGREE;
    FileName fileName, covName;

    readConfig(config, "inputfileTimeSplinesGravityfield", fileName,  Config::MUSTSET,  "{groopsDataDir}/", "");
    readConfig(config, "inputfileTimeSplinesCovariance",   covName,   Config::OPTIONAL, "",    "");
    readConfig(config, "minDegree",                        minDegree, Config::DEFAULT,  "0",   "");
    readConfig(config, "maxDegree",                        maxDegree, Config::OPTIONAL, "",    "");
    readConfig(config, "factor",                           factor,    Config::DEFAULT,  "1.0", "the result is multplied by this factor, set -1 to substract the field");
    if(isCreateSchema(config)) return;

    splinesFile.open(fileName, maxDegree, minDegree);
    hasCovariance = !covName.empty();
    if(hasCovariance)
      covarianceFile.open(covName, maxDegree, minDegree);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double GravityfieldTimeSplines::potential(const Time &time, const Vector3d &point) const
{
  return splinesFile.sphericalHarmonics(time, factor).potential(point);
}

/***********************************************/

Double GravityfieldTimeSplines::radialGradient(const Time &time, const Vector3d &point) const
{
  return splinesFile.sphericalHarmonics(time, factor).radialGradient(point);
}

/***********************************************/

Double GravityfieldTimeSplines::field(const Time &time, const Vector3d &point, const Kernel &kernel) const
{
  SphericalHarmonics harm = splinesFile.sphericalHarmonics(time, factor);
  // Convolution with kernel
  return inner(kernel.inverseCoefficients(point, harm.maxDegree()), harm.Yn(point, harm.maxDegree()));
}

/***********************************************/

Vector3d GravityfieldTimeSplines::gravity(const Time &time, const Vector3d &point) const
{
  return splinesFile.sphericalHarmonics(time, factor).gravity(point);
}

/***********************************************/

Tensor3d GravityfieldTimeSplines::gravityGradient(const Time &time, const Vector3d &point) const
{
  return splinesFile.sphericalHarmonics(time, factor).gravityGradient(point);
}

/***********************************************/

Vector3d GravityfieldTimeSplines::deformation(const Time &time, const Vector3d &point, Double gravity, const Vector &hn, const Vector &ln) const
{
  return splinesFile.sphericalHarmonics(time, factor).deformation(point, gravity, hn, ln);
}

/***********************************************/

void GravityfieldTimeSplines::deformation(const std::vector<Time> &time, const std::vector<Vector3d> &point, const std::vector<Double> &gravity,
                                          const Vector &hn, const Vector &ln, std::vector< std::vector<Vector3d> > &disp) const
{
  if((time.size()==0) || (point.size()==0))
    return;

  Matrix A;
  for(UInt i=0; i<time.size(); i++)
  {
    SphericalHarmonics harm = splinesFile.sphericalHarmonics(time.at(i), factor);
    Vector anm = harm.x();

    if(A.columns() < anm.rows())
      A = deformationMatrix(point, gravity, hn, ln, harm.GM(), harm.R(), harm.maxDegree());

    Vector x = A.column(0, anm.rows())*anm;
    for(UInt k=0; k<point.size(); k++)
    {
      disp.at(k).at(i).x() += x(3*k+0);
      disp.at(k).at(i).y() += x(3*k+1);
      disp.at(k).at(i).z() += x(3*k+2);
    }
  }
}

/***********************************************/

SphericalHarmonics GravityfieldTimeSplines::sphericalHarmonics(const Time &time, UInt maxDegree, UInt minDegree, Double GM, Double R) const
{
  try
  {
    if(time==Time())
      return SphericalHarmonics();
    SphericalHarmonics harm = splinesFile.sphericalHarmonics(time, factor).get(maxDegree, std::max(this->minDegree,minDegree), GM, R);
    if(!hasCovariance)
      return harm;

    Matrix C = covarianceFile.covariance(time, factor, maxDegree, std::max(this->minDegree,minDegree), harm.GM(), harm.R());
    if(C.size() == 0)
      return harm;
    Vector sigma2;
    if(C.getType() == Matrix::SYMMETRIC)
    {
      sigma2 = Vector(C.rows());
      for(UInt i=0; i<C.rows(); i++)
        sigma2(i) = C(i,i);
    }
    else
      sigma2 = C;

    Matrix sigma2cnm(harm.maxDegree()+1, Matrix::TRIANGULAR, Matrix::LOWER);
    Matrix sigma2snm(harm.maxDegree()+1, Matrix::TRIANGULAR, Matrix::LOWER);
    UInt idx = 0;
    for(UInt n=0; n<=harm.maxDegree(); n++)
    {
      sigma2cnm(n,0) = sigma2(idx++);
      for(UInt m=1; m<=n; m++)
      {
        sigma2cnm(n,m) = sigma2(idx++);
        sigma2snm(n,m) = sigma2(idx++);
      }
    }

    return SphericalHarmonics(harm.GM(), harm.R(), harm.cnm(), harm.snm(), sigma2cnm, sigma2snm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix GravityfieldTimeSplines::sphericalHarmonicsCovariance(const Time &time, UInt maxDegree, UInt minDegree, Double GM, Double R) const
{
  try
  {
    if((!hasCovariance) || (time==Time()))
      return Matrix();

    return covarianceFile.covariance(time, factor, maxDegree, std::max(this->minDegree,minDegree), GM, R);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GravityfieldTimeSplines::variance(const Time &time, const std::vector<Vector3d> &point, const Kernel &kernel, Matrix &D) const
{
  try
  {
    if((!hasCovariance) || (time==Time()))
      return;

    Double GM        = covarianceFile.GM();
    Double R         = covarianceFile.R();
    Matrix C         = covarianceFile.covariance(time, 1.);
    UInt   maxDegree = covarianceFile.maxDegree();

    // A = linear function from spherical harmonics to point values
    Matrix A(point.size(), C.rows());
    for(UInt k=0; k<point.size(); k++)
    {
      Matrix Cnm, Snm;
      SphericalHarmonics::CnmSnm(1/R * point.at(k), maxDegree, Cnm, Snm);
      Vector kn = kernel.inverseCoefficients(point.at(k), maxDegree);
      UInt idx = 0;
      for(UInt n=0; idx<A.columns(); n++)
      {
        A(k,idx++) = factor * kn(n) * GM/R * Cnm(n,0);
        for(UInt m=1; m<=n; m++)
        {
          A(k,idx++) = factor * kn(n) * GM/R * Cnm(n,m);
          A(k,idx++) = factor * kn(n) * GM/R * Snm(n,m);
        }
      }
    }

    if(C.getType() == Matrix::SYMMETRIC) // full covariance matrix?
      D += A * C * A.trans();
    else // only variances
    {
      for(UInt i=0; i<A.columns(); i++)
        A.column(i) *= sqrt(C(i,0));
      rankKUpdate(1., A.trans(), D);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
