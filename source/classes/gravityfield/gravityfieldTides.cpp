/***********************************************/
/**
* @file gravityfieldTides.cpp
*
* @brief Tidal forces computed with class Tide.
* @see Gravityfield
* @see Tides
*
* @author Torsten Mayer-Guerr
* @date 2004-11-01
*
*/
/***********************************************/

#include "base/import.h"
#include "base/sphericalHarmonics.h"
#include "config/config.h"
#include "classes/kernel/kernel.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/tides/tides.h"
#include "classes/gravityfield/gravityfield.h"
#include "classes/gravityfield/gravityfieldTides.h"

/***********************************************/

GravityfieldTides::GravityfieldTides(Config &config)
{
  try
  {
    readConfig(config, "tides",         tides,         Config::MUSTSET,  "",    "");
    readConfig(config, "earthRotation", earthRotation, Config::MUSTSET,  "",    "");
    readConfig(config, "ephemerides",   ephemerides,   Config::OPTIONAL, "jpl", "");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double GravityfieldTides::potential(const Time &time, const Vector3d &point) const
{
  return tides->potential(time, point, earthRotation->rotaryMatrix(time), earthRotation, ephemerides);
}

/***********************************************/

Double GravityfieldTides::radialGradient(const Time &time, const Vector3d &point) const
{
  return tides->radialGradient(time, point, earthRotation->rotaryMatrix(time), earthRotation, ephemerides);
}

/***********************************************/

Vector3d GravityfieldTides::gravity(const Time &time, const Vector3d &point) const
{
  return tides->acceleration(time, point, earthRotation->rotaryMatrix(time), earthRotation, ephemerides);
}

/***********************************************/

Tensor3d GravityfieldTides::gravityGradient(const Time &time, const Vector3d &point) const
{
  return tides->gradient(time, point, earthRotation->rotaryMatrix(time), earthRotation, ephemerides);
}

/***********************************************/

Vector3d GravityfieldTides::deformation(const Time &time, const Vector3d &point, Double gravity, const Vector &hn, const Vector &ln) const
{
  return tides->deformation(time, point, earthRotation->rotaryMatrix(time), earthRotation, ephemerides, gravity, hn, ln);
}

/***********************************************/

void GravityfieldTides::deformation(const std::vector<Time> &time, const std::vector<Vector3d> &point, const std::vector<Double> &gravity,
                                    const Vector &hn, const Vector &ln, std::vector<std::vector<Vector3d>> &disp) const
{
  std::vector<Rotary3d> rotEarth(time.size());
  for(UInt i=0; i<time.size(); i++)
    rotEarth[i]= earthRotation->rotaryMatrix(time[i]);
  tides->deformation(time, point, rotEarth, earthRotation, ephemerides, gravity, hn, ln, disp);
}

/***********************************************/

void GravityfieldTides::variance(const Time &/*time*/, const std::vector<Vector3d> &/*point*/, const Kernel &/*kernel*/, Matrix &/*D*/) const
{
  try
  {
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}


/***********************************************/

SphericalHarmonics GravityfieldTides::sphericalHarmonics(const Time &time, UInt maxDegree, UInt minDegree, Double GM, Double R) const
{
  try
  {
    return tides->sphericalHarmonics(time, earthRotation->rotaryMatrix(time), earthRotation, ephemerides, maxDegree, minDegree, GM, R);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
