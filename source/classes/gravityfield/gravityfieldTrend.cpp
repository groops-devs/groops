/***********************************************/
/**
* @file gravityfieldTrend.cpp
*
* @brief Gravityfield as trend.
* @see Gravityfield
*
* @author Torsten Mayer-Guerr
* @date 2007-06-15
*
*/
/***********************************************/

#include "base/import.h"
#include "base/sphericalHarmonics.h"
#include "config/config.h"
#include "classes/kernel/kernel.h"
#include "classes/gravityfield/gravityfield.h"
#include "classes/gravityfield/gravityfieldTrend.h"

/***********************************************/

GravityfieldTrend::GravityfieldTrend(Config &config)
{
  try
  {
    readConfig(config, "gravityfield", gravityfield, Config::MUSTSET, "",           "this field is multiplicated by (time-time0)/timeStep");
    readConfig(config, "timeStart",    time0,        Config::MUSTSET, STRING_J2000, "reference time");
    readConfig(config, "timeStep",     timeStep,     Config::MUSTSET, "365.25",     "");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double GravityfieldTrend::factor(Time time) const
{
  if(time == Time())
    return 0;
  return (time-time0).mjd()/timeStep.mjd();
}

/***********************************************/

Double GravityfieldTrend::potential(const Time &time, const Vector3d &point) const
{
  return factor(time) * gravityfield->potential(time, point);
}

/***********************************************/

Double GravityfieldTrend::radialGradient(const Time &time, const Vector3d &point) const
{
  return factor(time) * gravityfield->radialGradient(time, point);
}

/***********************************************/

Double GravityfieldTrend::field(const Time &time, const Vector3d &point, const Kernel &kernel) const
{
  return factor(time) * gravityfield->field(time, point, kernel);
}

/***********************************************/

Vector3d GravityfieldTrend::gravity(const Time &time, const Vector3d &point) const
{
  return factor(time) * gravityfield->gravity(time, point);
}

/***********************************************/

Tensor3d GravityfieldTrend::gravityGradient(const Time &time, const Vector3d &point) const
{
  return factor(time) * gravityfield->gravityGradient(time, point);
}

/***********************************************/

Vector3d GravityfieldTrend::deformation(const Time &time, const Vector3d &point, Double gravity, const Vector &hn, const Vector &ln) const
{
  return factor(time) * gravityfield->deformation(time, point, gravity, hn, ln);
}

/***********************************************/

void GravityfieldTrend::deformation(const std::vector<Time> &time, const std::vector<Vector3d> &point, const std::vector<Double> &gravity,
                                    const Vector &hn, const Vector &ln, std::vector<std::vector<Vector3d>> &disp) const
{
  if((time.size()==0) || (point.size()==0))
    return;

  std::vector<std::vector<Vector3d>> disp2(point.size());
  for(UInt k=0; k<point.size(); k++)
    disp2.at(k).resize(time.size());
  gravityfield->deformation(time, point, gravity, hn, ln, disp2);

  for(UInt i=0; i<time.size(); i++)
  {
    const Double f = factor(time.at(i));
    for(UInt k=0; k<point.size(); k++)
      disp.at(k).at(i) += f * disp2.at(k).at(i);
  }
}

/***********************************************/

void GravityfieldTrend::variance(const Time &/*time*/, const std::vector<Vector3d> &/*point*/, const Kernel &/*kernel*/, Matrix &/*D*/) const
{
}

/***********************************************/

SphericalHarmonics GravityfieldTrend::sphericalHarmonics(const Time &time, UInt maxDegree, UInt minDegree, Double GM, Double R) const
{
  return factor(time) * gravityfield->sphericalHarmonics(time, maxDegree, minDegree, GM, R);
}

/***********************************************/
