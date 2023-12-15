/***********************************************/
/**
* @file gravityfieldFromParametrization.cpp
*
* @brief Gravity field from parameter vector.
* @see Gravityfield
* @see ParametrizationGravity
*
* @author Torsten Mayer-Guerr
* @date 2003-10-31
*
*/
/***********************************************/

#include "base/import.h"
#include "base/sphericalHarmonics.h"
#include "config/config.h"
#include "files/fileMatrix.h"
#include "classes/kernel/kernel.h"
#include "classes/parametrizationGravity/parametrizationGravity.h"
#include "classes/gravityfield/gravityfield.h"
#include "classes/gravityfield/gravityfieldFromParametrization.h"

/***********************************************/

GravityfieldFromParametrization::GravityfieldFromParametrization(Config &config)
{
  try
  {
    UInt     rightSide, indexStart;
    FileName xName, WName, sigmaxName;

    renameDeprecatedConfig(config, "representation", "parametrization", date2time(2020, 6, 3));

    readConfig(config, "parametrization",   parametrization, Config::MUSTSET,  "",    "");
    readConfig(config, "inputfileSolution", xName,           Config::MUSTSET,  "",    "solution vector");
    readConfig(config, "inputfileSigmax",   sigmaxName,      Config::OPTIONAL, "",    "standards deviations or covariance matrix of the solution");
    readConfig(config, "indexStart",        indexStart,      Config::DEFAULT,  "0",   "position in the solution vector");
    readConfig(config, "rightSide",         rightSide,       Config::DEFAULT,  "0",   "if solution contains several right hand sides, select one");
    readConfig(config, "factor",            factor,          Config::DEFAULT,  "1.0", "the result is multiplied by this factor, set -1 to subtract the field");
    if(isCreateSchema(config)) return;

    Matrix mx;
    readFileMatrix(xName, mx);
    x = mx.slice(indexStart, rightSide, parametrization->parameterCount(), 1);

    if(!sigmaxName.empty())
    {
      Matrix A;
      readFileMatrix(sigmaxName, A);
      if(A.getType() == Matrix::SYMMETRIC)
      {
        // full covariance matrix
        C = A.slice(indexStart, indexStart, parametrization->parameterCount(), parametrization->parameterCount());
        sigma2x = Vector(C.rows());
        for(UInt i=0; i<sigma2x.rows(); i++)
          sigma2x(i) = C(i,i);
      }
      else
      {
        // only diagonal matrix
        sigma2x = A.row(indexStart, parametrization->parameterCount());
        for(UInt i=0; i<sigma2x.rows(); i++)
          sigma2x(i) *= sigma2x(i);
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double GravityfieldFromParametrization::potential(const Time &time, const Vector3d &point) const
{
  Matrix A(1, parametrization->parameterCount());
  parametrization->potential(time, point, A);
  return factor * (A*x)(0,0);
}

/***********************************************/

Double GravityfieldFromParametrization::radialGradient(const Time &time, const Vector3d &point) const
{
  Matrix A(1, parametrization->parameterCount());
  parametrization->radialGradient(time, point, A);
  return factor * (A*x)(0,0);
}

/***********************************************/

Double GravityfieldFromParametrization::field(const Time &time, const Vector3d &point, const Kernel &kernel) const
{
  Matrix A(1, parametrization->parameterCount());
  parametrization->field(time, point, kernel, A);
  return factor * (A*x)(0,0);
}

/***********************************************/

Vector3d GravityfieldFromParametrization::gravity(const Time &time, const Vector3d &point) const
{
  Matrix A(3, parametrization->parameterCount());
  parametrization->gravity(time, point, A);
  Vector l = A*x;
  Vector3d field;
  field.x() = l(0);
  field.y() = l(1);
  field.z() = l(2);
  return factor * field;
}

/***********************************************/

Tensor3d GravityfieldFromParametrization::gravityGradient(const Time &time, const Vector3d &point) const
{
  Matrix A(6, parametrization->parameterCount());
  parametrization->gravityGradient(time, point, A);
  Vector l = A*x;
  Tensor3d field;
  field.xx() = l(0);
  field.xy() = l(1);
  field.xz() = l(2);
  field.yy() = l(3);
  field.yz() = l(4);
  field.zz() = l(5);
  return factor * field;
}

/***********************************************/

Vector3d GravityfieldFromParametrization::deformation(const Time &time, const Vector3d &point, Double gravity, const Vector &hn, const Vector &ln) const
{
  Matrix A(3, parametrization->parameterCount());
  parametrization->deformation(time, point, gravity, hn, ln, A);
  const Vector l = A*x;
  return factor * Vector3d(l(0), l(1), l(2));
}

/***********************************************/

void GravityfieldFromParametrization::deformation(const std::vector<Time> &time, const std::vector<Vector3d> &point, const std::vector<Double> &gravity,
                                                 const Vector &hn, const Vector &ln, std::vector<std::vector<Vector3d>> &disp) const
{
  try
  {
    for(UInt i=0; i<time.size(); i++)
      for(UInt k=0; k<point.size(); k++)
        disp.at(k).at(i) += deformation(time.at(i), point.at(k), gravity.at(k), hn, ln);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GravityfieldFromParametrization::variance(const Time &time, const std::vector<Vector3d> &point, const Kernel &kernel, Matrix &D) const
{
  try
  {
    Matrix A(point.size(), parametrization->parameterCount());
    for(UInt i=0; i<point.size(); i++)
      parametrization->field(time, point.at(i), kernel, A.row(i));
    if(C.size())
      D += A*C*A.trans();
    else
    {
      for(UInt i=0; i<A.columns(); i++)
        A.column(i) *= sqrt(sigma2x(i));
      rankKUpdate(1., A.trans(), D);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}


/***********************************************/

SphericalHarmonics GravityfieldFromParametrization::sphericalHarmonics(const Time &time, UInt maxDegree, UInt minDegree, Double GM, Double R) const
{
  try
  {
    if(sigma2x.rows()!=0)
      return factor * parametrization->sphericalHarmonics(time, x, sigma2x, maxDegree).get(maxDegree, minDegree, GM, R);
    else
      return factor * parametrization->sphericalHarmonics(time, x, maxDegree).get(maxDegree, minDegree, GM, R);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
