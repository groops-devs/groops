/***********************************************/
/**
* @file parametrizationGravityRadialBasis.cpp
*
* @brief Radial basis functions.
* @see ParametrizationGravity
*
* @author Torsten Mayer-Guerr
* @date 2003-03-10
*
*/
/***********************************************/

#include "base/import.h"
#include "base/sphericalHarmonics.h"
#include "config/config.h"
#include "classes/grid/grid.h"
#include "classes/kernel/kernel.h"
#include "classes/parametrizationGravity/parametrizationGravity.h"
#include "classes/parametrizationGravity/parametrizationGravityRadialBasis.h"

/***********************************************/

ParametrizationGravityRadialBasis::ParametrizationGravityRadialBasis(Config &config)
{
  try
  {
    GridPtr grid;

    readConfig(config, "kernel",kernel, Config::MUSTSET, "", "shape of the radial basis function");
    readConfig(config, "grid",  grid,   Config::MUSTSET, "", "nodal point distribution");
    if(isCreateSchema(config)) return;

    sourcePoint = grid->points();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ParametrizationGravityRadialBasis::parameterName(std::vector<ParameterName> &name) const
{
  for(UInt i=0; i<sourcePoint.size(); i++)
    name.push_back(ParameterName("", "radialBasis."+i%"%i"s+"."+sourcePoint.size()%"%i"s));
}

/***********************************************/

void ParametrizationGravityRadialBasis::field(const Time &/*time*/, const Vector3d &point, const Kernel &kernel2, MatrixSliceRef A) const
{
  for(UInt i=0; i<sourcePoint.size(); i++)
    A(0,i) = kernel2.inverseKernel(point, sourcePoint.at(i), *kernel);
}

/***********************************************/

void ParametrizationGravityRadialBasis::potential(const Time &/*time*/, const Vector3d &point, MatrixSliceRef A) const
{
  for(UInt i=0; i<sourcePoint.size(); i++)
    A(0,i) = kernel->kernel(point, sourcePoint.at(i));
}

/***********************************************/

void ParametrizationGravityRadialBasis::radialGradient(const Time &/*time*/, const Vector3d &point, MatrixSliceRef A) const
{
  for(UInt i=0; i<sourcePoint.size(); i++)
    A(0,i) = kernel->radialDerivative(point, sourcePoint.at(i));
}

/***********************************************/

void ParametrizationGravityRadialBasis::gravity(const Time &/*time*/, const Vector3d &point, MatrixSliceRef A) const
{
  for(UInt i=0; i<sourcePoint.size(); i++)
  {
    Vector3d field = kernel->gradient(point, sourcePoint.at(i));
    A(0,i) = field.x();
    A(1,i) = field.y();
    A(2,i) = field.z();
  }
}

/***********************************************/

void ParametrizationGravityRadialBasis::gravityGradient(const Time &/*time*/, const Vector3d &point, MatrixSliceRef A) const
{
  for(UInt i=0; i<sourcePoint.size(); i++)
  {
    Tensor3d field = kernel->gradientGradient(point, sourcePoint.at(i));
    A(0,i) = field.xx();
    A(1,i) = field.xy();
    A(2,i) = field.xz();
    A(3,i) = field.yy();
    A(4,i) = field.yz();
    A(5,i) = field.zz();
  }
}

/***********************************************/

void ParametrizationGravityRadialBasis::deformation(const Time &/*time*/, const Vector3d &/*point*/, Double /*gravity*/, const Vector &/*hn*/, const Vector &/*ln*/, MatrixSliceRef /*A*/) const
{
  try
  {
    throw(Exception("not implemented"));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

SphericalHarmonics ParametrizationGravityRadialBasis::sphericalHarmonics(const Time &/*time*/, const Vector &x, UInt maxDegree) const
{
  try
  {
    // Determine max. degree of the kernel
    if(maxDegree == INFINITYDEGREE)
      maxDegree = kernel->maxDegree();
    if(maxDegree == INFINITYDEGREE)
      throw(Exception("In ParametrizationGravityRadialBasis::sphericalHarmonics: INFINITYDEGREE requested"));

    Matrix cnm(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    Matrix snm(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);

    for(UInt i=0; i<sourcePoint.size(); i++)
    {
      Vector coeff = kernel->coefficients(sourcePoint.at(i), maxDegree);

      Double Ri = sourcePoint.at(i).r();
      Matrix Cnm, Snm;
      SphericalHarmonics::CnmSnm((1/Ri) * sourcePoint.at(i), maxDegree, Cnm, Snm);

      Double factor = DEFAULT_R/DEFAULT_GM * Ri/DEFAULT_R;
      for(UInt n=0; n<=maxDegree; n++)
      {
        for(UInt m=0; m<=n; m++)
        {
          cnm(n,m) += factor * coeff(n) * Cnm(n,m) * x(i);
          snm(n,m) += factor * coeff(n) * Snm(n,m) * x(i);
        }
        factor *= Ri/DEFAULT_R;
      }
    }

    return SphericalHarmonics(DEFAULT_GM, DEFAULT_R, cnm, snm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

SphericalHarmonics ParametrizationGravityRadialBasis::sphericalHarmonics(const Time &/*time*/, const Vector &x, const Vector &sigma2x, UInt maxDegree) const
{
  try
  {
    // Determine max. degree of the kernel
    if(maxDegree == INFINITYDEGREE)
      maxDegree = kernel->maxDegree();
    if(maxDegree == INFINITYDEGREE)
      throw(Exception("In ParametrizationGravityRadialBasis::sphericalHarmonics: INFINITYDEGREE requested"));

    Matrix cnm      (maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    Matrix snm      (maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    Matrix sigma2cnm(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    Matrix sigma2snm(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);

    // mit Fehlerfortplanzung
    for(UInt i=0; i<sourcePoint.size(); i++)
    {
      Vector coeff = kernel->coefficients(sourcePoint.at(i), maxDegree);

      Double Ri = sourcePoint.at(i).r();
      Matrix Cnm, Snm;
      SphericalHarmonics::CnmSnm((1/Ri) * sourcePoint.at(i), maxDegree, Cnm, Snm);

      Double factor = DEFAULT_R/DEFAULT_GM * Ri/DEFAULT_R;
      for(UInt n=0; n<=maxDegree; n++)
      {
        for(UInt m=0; m<=n; m++)
        {
          cnm(n,m)       += factor * coeff(n) * Cnm(n,m) * x(i);
          snm(n,m)       += factor * coeff(n) * Snm(n,m) * x(i);
          sigma2cnm(n,m) += pow(factor * coeff(n) * Cnm(n,m), 2) * sigma2x(i);
          sigma2snm(n,m) += pow(factor * coeff(n) * Snm(n,m), 2) * sigma2x(i);
        }
        factor *= Ri/DEFAULT_R;
      }
    }

    return SphericalHarmonics(DEFAULT_GM, DEFAULT_R, cnm, snm, sigma2cnm, sigma2snm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
