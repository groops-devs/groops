/***********************************************/
/**
* @file parametrizationGravityTemporal.cpp
*
* @brief Time variable gravity field.
* @see ParametrizationGravity
*
* @author Torsten Mayer-Guerr
* @date 2011-12-28
*
*/
/***********************************************/

#include "base/import.h"
#include "base/sphericalHarmonics.h"
#include "config/config.h"
#include "classes/kernel/kernel.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "classes/parametrizationGravity/parametrizationGravity.h"
#include "parametrizationGravityTemporal.h"

/***********************************************/

ParametrizationGravityTemporal::ParametrizationGravityTemporal(Config &config)
{
  try
  {
    renameDeprecatedConfig(config, "representation", "parametrizationGravity",  date2time(2020, 6, 3));
    renameDeprecatedConfig(config, "temporal",       "parametrizationTemporal", date2time(2020, 6, 3));

    readConfig(config, "parametrizationGravity",  spatial,  Config::MUSTSET, "", "");
    readConfig(config, "parametrizationTemporal", temporal, Config::MUSTSET, "", "");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ParametrizationGravityTemporal::setInterval(const Time &timeStart, const Time &timeEnd)
{
  try
  {
    temporal->setInterval(timeStart, timeEnd);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ParametrizationGravityTemporal::parameterName(std::vector<ParameterName> &name) const
{
  try
  {
    std::vector<ParameterName> baseName;
    spatial->parameterName(baseName);
    temporal->parameterName(baseName, name);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ParametrizationGravityTemporal::field(const Time &time, const Vector3d &point, const Kernel &kernel2, MatrixSliceRef A) const
{
  try
  {
    Matrix B(1,spatial->parameterCount());
    spatial->field(time, point, kernel2, B);
    temporal->designMatrix(time, B, A);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ParametrizationGravityTemporal::potential(const Time &time, const Vector3d &point, MatrixSliceRef A) const
{
  try
  {
    Matrix B(1,spatial->parameterCount());
    spatial->potential(time, point, B);
    temporal->designMatrix(time, B, A);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ParametrizationGravityTemporal::radialGradient(const Time &time, const Vector3d &point, MatrixSliceRef A) const
{
  try
  {
    Matrix B(1,spatial->parameterCount());
    spatial->radialGradient(time, point, B);
    temporal->designMatrix(time, B, A);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ParametrizationGravityTemporal::gravity(const Time &time, const Vector3d &point, MatrixSliceRef A) const
{
  try
  {
    Matrix B(3,spatial->parameterCount());
    spatial->gravity(time, point, B);
    temporal->designMatrix(time, B, A);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ParametrizationGravityTemporal::gravityGradient(const Time &time, const Vector3d &point, MatrixSliceRef A) const
{
  try
  {
    Matrix B(6,spatial->parameterCount());
    spatial->gravityGradient(time, point, B);
    temporal->designMatrix(time, B, A);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ParametrizationGravityTemporal::deformation(const Time &time, const Vector3d &point, Double gravity, const Vector &hn, const Vector &ln, MatrixSliceRef A) const
{
  try
  {
    Matrix B(3,spatial->parameterCount());
    spatial->deformation(time, point, gravity, hn, ln, A);
    temporal->designMatrix(time, B, A);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

SphericalHarmonics ParametrizationGravityTemporal::sphericalHarmonics(const Time &time, const Vector &x, UInt maxDegree) const
{
  try
  {
    if(time==Time())
      return ::SphericalHarmonics();

    std::vector<UInt>   index;
    std::vector<Double> factor;
    temporal->factors(time, index, factor);

    SphericalHarmonics harm;
    UInt count = spatial->parameterCount();
    for(UInt i=0; i<index.size(); i++)
      harm  += factor.at(i) * spatial->sphericalHarmonics(time, x.row(index.at(i)*count,count), maxDegree);

    return harm;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

SphericalHarmonics ParametrizationGravityTemporal::sphericalHarmonics(const Time &time, const Vector &x, const Vector &sigma2x, UInt maxDegree) const
{
  try
  {
    if(time==Time())
      return ::SphericalHarmonics();

    std::vector<UInt>   index;
    std::vector<Double> factor;
    temporal->factors(time, index, factor);

    SphericalHarmonics harm;
    UInt count = spatial->parameterCount();
    for(UInt i=0; i<index.size(); i++)
      harm  += factor.at(i) * spatial->sphericalHarmonics(time, x.row(index.at(i)*count,count), sigma2x.row(index.at(i)*count,count), maxDegree);

    return harm;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
