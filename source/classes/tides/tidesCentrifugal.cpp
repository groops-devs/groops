/***********************************************/
/**
* @file tidesCentrifugal.cpp
*
* @brief Centrifugal force.
* Actual centrifugal force computed from Earth rotation at given time.
* @see Tides
*
* @author Torsten Mayer-Guerr
* @date 2011-10-19
*
*/
/***********************************************/

#include "base/import.h"
#include "base/sphericalHarmonics.h"
#include "config/config.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/tides/tides.h"
#include "classes/tides/tidesCentrifugal.h"

/***********************************************/

TidesCentrifugal::TidesCentrifugal(Config &config)
{
  try
  {
    readConfig(config, "factor", factor, Config::DEFAULT,  "1.0", "the result is multiplied by this factor, set -1 to subtract the field");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double TidesCentrifugal::potential(const Time &time, const Vector3d &point, const Rotary3d &/*rotEarth*/, EarthRotationPtr rotation, EphemeridesPtr /*ephemerides*/) const
{
  Vector3d Omega = rotation->rotaryAxis((time==Time()) ? mjd2time(J2000) : time);
  return factor * 0.5 * crossProduct(Omega, point).quadsum();
}

/***********************************************/

Double TidesCentrifugal::radialGradient(const Time &time, const Vector3d &point, const Rotary3d &rotEarth, EarthRotationPtr rotation, EphemeridesPtr ephemerides) const
{
  // df/dr = df/dpoint * dpoint/dr
  // with dpoint/dr = (1/norm(point))*point
  return inner(gravity(time, point, rotEarth, rotation, ephemerides), point)/point.r();
}

/***********************************************/

Vector3d TidesCentrifugal::gravity(const Time &time, const Vector3d &point, const Rotary3d &/*rotEarth*/, EarthRotationPtr rotation, EphemeridesPtr /*ephemerides*/) const
{
  Vector3d Omega = rotation->rotaryAxis((time==Time()) ? mjd2time(J2000) : time);
  return -factor * crossProduct(Omega, crossProduct(Omega, point));
}

/***********************************************/

Tensor3d TidesCentrifugal::gravityGradient(const Time &time, const Vector3d &/*point*/, const Rotary3d &/*rotEarth*/, EarthRotationPtr rotation, EphemeridesPtr /*ephemerides*/) const
{
  Vector3d Omega = rotation->rotaryAxis((time==Time()) ? mjd2time(J2000) : time);
  Tensor3d T;
  T.xx() = pow(Omega.y(),2)+pow(Omega.z(),2);
  T.yy() = pow(Omega.x(),2)+pow(Omega.z(),2);
  T.zz() = pow(Omega.x(),2)+pow(Omega.y(),2);
  T.xy() = -Omega.x() * Omega.y();
  T.xz() = -Omega.x() * Omega.z();
  T.yz() = -Omega.y() * Omega.z();
  return factor * T;
}

/***********************************************/

Vector3d TidesCentrifugal::deformation(const Time &/*time*/, const Vector3d &/*point*/, const Rotary3d &/*rotEarth*/, EarthRotationPtr /*earthRotation*/, EphemeridesPtr /*ephemerides*/,
                                       Double /*gravity*/, const Vector &/*hn*/, const Vector &/*ln*/) const
{
  return Vector3d(0,0,0);
}

/***********************************************/

void TidesCentrifugal::deformation(const std::vector<Time> &/*time*/, const std::vector<Vector3d> &/*point*/, const std::vector<Rotary3d> &/*rotEarth*/,
                                   EarthRotationPtr /*rotation*/, EphemeridesPtr /*ephemerides*/, const std::vector<Double> &/*gravity*/, const Vector &/*hn*/, const Vector &/*ln*/,
                                   std::vector<std::vector<Vector3d>> &/*disp*/) const
{
}

/***********************************************/

SphericalHarmonics TidesCentrifugal::sphericalHarmonics(const Time &/*time*/, const Rotary3d &/*rotEarth*/, EarthRotationPtr /*earthRotation*/, EphemeridesPtr /*ephemerides*/,
                                                        UInt /*maxDegree*/, UInt /*minDegree*/, Double /*GM*/, Double /*R*/) const
{
  return SphericalHarmonics();
}

/***********************************************/
