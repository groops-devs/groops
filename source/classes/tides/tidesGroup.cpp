/***********************************************/
/**
* @file tidesGroup.cpp
*
* @brief Group.
* @see Tides
*
* @author Torsten Mayer-Guerr
* @date 2023-11-13
*
*/
/***********************************************/

#include "base/import.h"
#include "base/sphericalHarmonics.h"
#include "config/config.h"
#include "files/fileMeanPolarMotion.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/tides/tides.h"
#include "classes/tides/tidesGroup.h"

/***********************************************/

TidesGroup::TidesGroup(Config &config)
{
  try
  {
    readConfig(config, "tides",  tides,  Config::DEFAULT, "", "");
    readConfig(config, "factor", factor, Config::DEFAULT, "1.0", "the result is multiplied by this factor");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}


/***********************************************/

Double TidesGroup::potential(const Time &time, const Vector3d &point, const Rotary3d &rotEarth, EarthRotationPtr rotation, EphemeridesPtr ephemerides) const
{
  try
  {
    return factor * tides->potential(time, point, rotEarth, rotation, ephemerides);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double TidesGroup::radialGradient(const Time &time, const Vector3d &point, const Rotary3d &rotEarth, EarthRotationPtr rotation, EphemeridesPtr ephemerides) const
{
  try
  {
    return factor * tides->radialGradient(time, point, rotEarth, rotation, ephemerides);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector3d TidesGroup::gravity(const Time &time, const Vector3d &point, const Rotary3d &rotEarth, EarthRotationPtr rotation, EphemeridesPtr ephemerides) const
{
  try
  {
    return factor * tides->acceleration(time, point, rotEarth, rotation, ephemerides);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Tensor3d TidesGroup::gravityGradient(const Time &time, const Vector3d &point, const Rotary3d &rotEarth, EarthRotationPtr rotation, EphemeridesPtr ephemerides) const
{
  try
  {
    return factor * tides->gradient(time, point, rotEarth, rotation, ephemerides);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector3d TidesGroup::deformation(const Time &time, const Vector3d &point, const Rotary3d &rotEarth, EarthRotationPtr rotation, EphemeridesPtr ephemerides,
                                 Double gravity, const Vector &hn, const Vector &ln) const
{
  try
  {
    return factor * tides->deformation(time, point, rotEarth, rotation, ephemerides, gravity, hn, ln);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void TidesGroup::deformation(const std::vector<Time> &time, const std::vector<Vector3d> &point, const std::vector<Rotary3d> &rotEarth,
                             EarthRotationPtr rotation, EphemeridesPtr ephemerides, const std::vector<Double> &gravity, const Vector &hn, const Vector &ln,
                             std::vector<std::vector<Vector3d>> &disp) const
{
  try
  {
    tides->deformation(time, point, rotEarth, rotation, ephemerides, gravity, hn, ln, disp);
    if(factor != 1.)
      for(auto &ds : disp)
        for(auto &d : ds)
          d *= factor;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

SphericalHarmonics TidesGroup::sphericalHarmonics(const Time &time, const Rotary3d &rotEarth, EarthRotationPtr rotation, EphemeridesPtr ephemerides,
                                                  UInt maxDegree, UInt minDegree, Double GM, Double R) const
{
  try
  {
    return factor * tides->sphericalHarmonics(time, rotEarth, rotation, ephemerides, maxDegree, minDegree, GM, R);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/
