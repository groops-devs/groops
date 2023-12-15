/***********************************************/
/**
* @file gravityfieldEarthquakeOscillation.cpp
*
* @brief Earthquake oscillation.
* @see Gravityfield
*
* @author Saniya Behzadpour
* @date 2017-06-07
*
*/
/***********************************************/

#include "base/import.h"
#include "base/sphericalHarmonics.h"
#include "config/config.h"
#include "files/fileMatrix.h"
#include "classes/kernel/kernel.h"
#include "classes/gravityfield/gravityfield.h"
#include "classes/gravityfield/gravityfieldEarthquakeOscillation.h"

/***********************************************/

GravityfieldEarthquakeOscillation::GravityfieldEarthquakeOscillation(Config &config)
{
  try
  {
    FileName xName;

    readConfig(config, "inputCoefficientMatrix", xName,     Config::MUSTSET,  "",  "oscillation model parameters");
    readConfig(config, "time0",                  time0,     Config::MUSTSET,  "",  "the time earthquake happened");
    readConfig(config, "minDegree",              minDegree, Config::DEFAULT,  "2", "");
    readConfig(config, "maxDegree",              maxDegree, Config::DEFAULT,  "2", "");
    readConfig(config, "GM",                     GM,        Config::DEFAULT,  STRING_DEFAULT_GM, "Geocentric gravitational constant");
    readConfig(config, "R",                      R,         Config::DEFAULT,  STRING_DEFAULT_R,  "reference radius");
    if(isCreateSchema(config)) return;

    Matrix mx;
    readFileMatrix(xName, mx);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}
/***********************************************/

SphericalHarmonics GravityfieldEarthquakeOscillation::timeVariableCoefficients(const Time &time) const
{
  try
  {
    if(time < time0)
      return SphericalHarmonics();

    const Double c = 2*PI*(time-time0).seconds();

    Matrix cnm(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    Matrix snm(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    for(UInt i=0; i<mx.rows(); i++)
    {
      cnm(mx(i,1), mx(i,2)) += mx(i,3) * (1-cos(c/mx(i,5)) * exp(-c/mx(i,5)/mx(i,6)/2));
      snm(mx(i,1), mx(i,2)) += mx(i,4) * (1-cos(c/mx(i,5)) * exp(-c/mx(i,5)/mx(i,6)/2));
    }

    return SphericalHarmonics(GM, R, cnm, snm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}
/***********************************************/

Double GravityfieldEarthquakeOscillation::potential(const Time &time, const Vector3d &point) const
{
  return timeVariableCoefficients(time).potential(point);
}

/***********************************************/

Double GravityfieldEarthquakeOscillation::radialGradient(const Time &time, const Vector3d &point) const
{
  return timeVariableCoefficients(time).radialGradient(point);
}

/***********************************************/

Double GravityfieldEarthquakeOscillation::field(const Time &time, const Vector3d &point, const Kernel &kernel) const
{
  SphericalHarmonics harmonics = timeVariableCoefficients(time);
  return inner(kernel.inverseCoefficients(point, harmonics.maxDegree(), harmonics.isInterior()), harmonics.Yn(point, harmonics.maxDegree()));
}

/***********************************************/

Vector3d GravityfieldEarthquakeOscillation::gravity(const Time &time, const Vector3d &point) const
{
  return timeVariableCoefficients(time).gravity(point);
}

/***********************************************/

Tensor3d GravityfieldEarthquakeOscillation::gravityGradient(const Time &time, const Vector3d &point) const
{
  return timeVariableCoefficients(time).gravityGradient(point);
}

/***********************************************/

Vector3d GravityfieldEarthquakeOscillation::deformation(const Time &time, const Vector3d &point, Double gravity, const Vector &hn, const Vector &ln) const
{
  return timeVariableCoefficients(time).deformation(point, gravity, hn, ln);
}

/***********************************************/

void GravityfieldEarthquakeOscillation::deformation(const std::vector<Time> &time, const std::vector<Vector3d> &point, const std::vector<Double> &gravity,
                                                    const Vector &hn, const Vector &ln, std::vector<std::vector<Vector3d>> &disp) const
{
  for(UInt i=0; i<time.size(); i++)
    for(UInt k=0; k<point.size(); k++)
      disp.at(k).at(i) += deformation(time.at(i), point.at(k), gravity.at(k), hn, ln);
}

/***********************************************/

SphericalHarmonics GravityfieldEarthquakeOscillation::sphericalHarmonics(const Time &time, UInt maxDegree, UInt minDegree, Double GM, Double R) const
{
  return timeVariableCoefficients(time).get(maxDegree, std::max(this->minDegree, minDegree), GM, R);
}

/***********************************************/

void GravityfieldEarthquakeOscillation::variance(const Time &/*time*/, const std::vector<Vector3d> &/*point*/, const Kernel &/*kernel*/, Matrix &/*D*/) const
{
}

/***********************************************/
