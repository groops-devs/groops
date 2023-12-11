/***********************************************/
/**
* @file gravityfieldOscillation.cpp
*
* @brief Gravityfield as oscillation.
* @see Gravityfield
*
* @author Torsten Mayer-Guerr
* @date 2015-06-07
*
*/
/***********************************************/

#include "base/import.h"
#include "base/sphericalHarmonics.h"
#include "config/config.h"
#include "classes/kernel/kernel.h"
#include "classes/gravityfield/gravityfield.h"
#include "classes/gravityfield/gravityfieldOscillation.h"

/***********************************************/

GravityfieldOscillation::GravityfieldOscillation(Config &config)
{
  try
  {
    readConfig(config, "gravityfieldCos", gravityfieldCos, Config::MUSTSET, "", "multiplicated by cos(2pi/T(time-time0))");
    readConfig(config, "gravityfieldSin", gravityfieldSin, Config::MUSTSET, "", "multiplicated by sin(2pi/T(time-time0))");
    readConfig(config, "time0",           time0,           Config::MUSTSET, STRING_J2000, "reference time");
    readConfig(config, "period",          timePeriod,      Config::MUSTSET, "365.25",     "[day]");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double GravityfieldOscillation::potential(const Time &time, const Vector3d &point) const
{
  if(time == Time())
    return 0;
  const Double omega = 2*PI/timePeriod.mjd()*(time-time0).mjd();
  return cos(omega) * gravityfieldCos->potential(time, point) + sin(omega) * gravityfieldSin->potential(time, point);
}

/***********************************************/

Double GravityfieldOscillation::radialGradient(const Time &time, const Vector3d &point) const
{
  if(time == Time())
    return 0;
  const Double omega = 2*PI/timePeriod.mjd()*(time-time0).mjd();
  return cos(omega) * gravityfieldCos->radialGradient(time, point) +
         sin(omega) * gravityfieldSin->radialGradient(time, point);
}

/***********************************************/

Double GravityfieldOscillation::field(const Time &time, const Vector3d &point, const Kernel &kernel) const
{
  if(time == Time())
    return 0;
  const Double omega = 2*PI/timePeriod.mjd()*(time-time0).mjd();
  return cos(omega) * gravityfieldCos->field(time, point, kernel) +
         sin(omega) * gravityfieldSin->field(time, point, kernel);
}

/***********************************************/

Vector3d GravityfieldOscillation::gravity(const Time &time, const Vector3d &point) const
{
  if(time == Time())
    return Vector3d();
  const Double omega = 2*PI/timePeriod.mjd()*(time-time0).mjd();
  return cos(omega) * gravityfieldCos->gravity(time, point) +
         sin(omega) * gravityfieldSin->gravity(time, point);
}

/***********************************************/

Tensor3d GravityfieldOscillation::gravityGradient(const Time &time, const Vector3d &point) const
{
  if(time == Time())
    return Tensor3d();
  const Double omega = 2*PI/timePeriod.mjd()*(time-time0).mjd();
  return cos(omega) * gravityfieldCos->gravityGradient(time, point) +
         sin(omega) * gravityfieldSin->gravityGradient(time, point);
}

/***********************************************/

Vector3d GravityfieldOscillation::deformation(const Time &time, const Vector3d &point, Double gravity, const Vector &hn, const Vector &ln) const
{
  if(time == Time())
    return Vector3d();
  const Double omega = 2*PI/timePeriod.mjd()*(time-time0).mjd();
  return cos(omega) * gravityfieldCos->deformation(time, point, gravity, hn, ln) +
         sin(omega) * gravityfieldSin->deformation(time, point, gravity, hn, ln);
}

/***********************************************/

void GravityfieldOscillation::deformation(const std::vector<Time> &time, const std::vector<Vector3d> &point, const std::vector<Double> &gravity,
                                    const Vector &hn, const Vector &ln, std::vector<std::vector<Vector3d>> &disp) const
{
  if((time.size()==0) || (point.size()==0))
    return;

  std::vector<std::vector<Vector3d>> dispCos(point.size());
  for(UInt k=0; k<point.size(); k++)
    dispCos.at(k).resize(time.size());
  gravityfieldCos->deformation(time, point, gravity, hn, ln, dispCos);

  std::vector<std::vector<Vector3d>> dispSin(point.size());
  for(UInt k=0; k<point.size(); k++)
    dispSin.at(k).resize(time.size());
  gravityfieldCos->deformation(time, point, gravity, hn, ln, dispSin);

  for(UInt i=0; i<time.size(); i++)
  {
    const Double omega = 2*PI/timePeriod.mjd()*(time.at(i)-time0).mjd();
    for(UInt k=0; k<point.size(); k++)
      disp.at(k).at(i) += cos(omega) * dispCos.at(k).at(i) + sin(omega) * dispSin.at(k).at(i);
  }
}

/***********************************************/

void GravityfieldOscillation::variance(const Time &/*time*/, const std::vector<Vector3d> &/*point*/, const Kernel &/*kernel*/, Matrix &/*D*/) const
{
}

/***********************************************/

SphericalHarmonics GravityfieldOscillation::sphericalHarmonics(const Time &time, UInt maxDegree, UInt minDegree, Double GM, Double R) const
{
  if(time == Time())
    return SphericalHarmonics();
  const Double omega = 2*PI/timePeriod.mjd()*(time-time0).mjd();
  return cos(omega) * gravityfieldCos->sphericalHarmonics(time, maxDegree, minDegree, GM, R) +
         sin(omega) * gravityfieldSin->sphericalHarmonics(time, maxDegree, minDegree, GM, R);
}

/***********************************************/
