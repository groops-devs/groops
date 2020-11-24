/***********************************************/
/**
* @file parametrizationGravityLinearTransformation.cpp
*
* @brief Gravity field parametrization based on the linear transformation of another parametrizationGravity.
* @see ParametrizationGravity
*
* @author Andreas Kvas
* @date 2020-06-28
*
*/
/***********************************************/

#include "base/import.h"
#include "base/sphericalHarmonics.h"
#include "config/config.h"
#include "classes/kernel/kernel.h"
#include "files/fileMatrix.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "classes/parametrizationGravity/parametrizationGravity.h"
#include "parametrizationGravityLinearTransformation.h"

/***********************************************/

ParametrizationGravityLinearTransformation::ParametrizationGravityLinearTransformation(Config &config)
{
  try
  {
    FileName fileNameTransformationMatrix;

    readConfig(config, "parametrizationGravitySource",  source, Config::MUSTSET, "", "");
    readConfig(config, "inputfileTransformationMatrix", fileNameTransformationMatrix, Config::MUSTSET, "", "transformation matrix from target to source parametrization (rows of this matrix must coincide with the parameter count of the source parametrization)");
    if(isCreateSchema(config)) return;

    readFileMatrix(fileNameTransformationMatrix, transformationMatrix);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ParametrizationGravityLinearTransformation::parameterName(std::vector<ParameterName> &name) const
{
  try
  {
    for(UInt i=0; i<parameterCount(); i++)
      name.push_back(ParameterName("", "transformedParameter."+i%"%i"s+"."+parameterCount()%"%i"s));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ParametrizationGravityLinearTransformation::field(const Time &time, const Vector3d &point, const Kernel &kernel2, MatrixSliceRef A) const
{
  try
  {
    Matrix B(1,source->parameterCount());
    source->field(time, point, kernel2, B);
    matMult(1.0, B, transformationMatrix, A);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ParametrizationGravityLinearTransformation::potential(const Time &time, const Vector3d &point, MatrixSliceRef A) const
{
  try
  {
    Matrix B(1,source->parameterCount());
    source->potential(time, point, B);
    matMult(1.0, B, transformationMatrix, A);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ParametrizationGravityLinearTransformation::radialGradient(const Time &time, const Vector3d &point, MatrixSliceRef A) const
{
  try
  {
    Matrix B(1,source->parameterCount());
    source->radialGradient(time, point, B);
    matMult(1.0, B, transformationMatrix, A);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ParametrizationGravityLinearTransformation::gravity(const Time &time, const Vector3d &point, MatrixSliceRef A) const
{
  try
  {
    Matrix B(3,source->parameterCount());
    source->gravity(time, point, B);
    matMult(1.0, B, transformationMatrix, A);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ParametrizationGravityLinearTransformation::gravityGradient(const Time &time, const Vector3d &point, MatrixSliceRef A) const
{
  try
  {
    Matrix B(6,source->parameterCount());
    source->gravityGradient(time, point, B);
    matMult(1.0, B, transformationMatrix, A);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ParametrizationGravityLinearTransformation::deformation(const Time &time, const Vector3d &point, Double gravity, const Vector &hn, const Vector &ln, MatrixSliceRef A) const
{
  try
  {
    Matrix B(3,source->parameterCount());
    source->deformation(time, point, gravity, hn, ln, A);
    matMult(1.0, B, transformationMatrix, A);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

SphericalHarmonics ParametrizationGravityLinearTransformation::sphericalHarmonics(const Time &time, const Vector &x, UInt maxDegree) const
{
  try
  {
    return source->sphericalHarmonics(time, transformationMatrix*x, maxDegree);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

SphericalHarmonics ParametrizationGravityLinearTransformation::sphericalHarmonics(const Time &time, const Vector &x, const Vector &sigma2x, UInt maxDegree) const
{
  try
  {
    // variance propagation
    Vector sigma2y(transformationMatrix.rows());
    for(UInt i=0; i<transformationMatrix.rows(); i++)
      for(UInt k=0; k<transformationMatrix.columns(); k++)
        sigma2y(i) = std::pow(transformationMatrix(i, k), 2) * sigma2x(k);

    return source->sphericalHarmonics(time, transformationMatrix*x, sigma2y, maxDegree);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
