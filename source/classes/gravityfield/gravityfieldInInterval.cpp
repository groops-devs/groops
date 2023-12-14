/***********************************************/
/**
* @file gravityfieldInInterval.cpp
*
* @brief Gravity fields which is not zero during a specific time span only.
* @see Gravityfield
*
* @author Torsten Mayer-Guerr
* @date 2017-11-06
*
*/
/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "inputOutput/logging.h"
#include "classes/gravityfield/gravityfield.h"
#include "classes/gravityfield/gravityfieldInInterval.h"

/***********************************************/

GravityfieldInInterval::GravityfieldInInterval(Config &config)
{
  try
  {
    readConfig(config, "gravityfield", gravityfield, Config::MUSTSET, "", "");
    readConfig(config, "timeStart",    timeStart,    Config::MUSTSET, "", "first point in time");
    readConfig(config, "timeEnd",      timeEnd,      Config::MUSTSET, "", "last point in time will be less or equal timeEnd");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double GravityfieldInInterval::potential(const Time &time, const Vector3d &point) const
{
  try
  {
    return time.isInInterval(timeStart, timeEnd) ? gravityfield->potential(time, point) : 0;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double GravityfieldInInterval::radialGradient(const Time &time, const Vector3d &point) const
{
  try
  {
    return time.isInInterval(timeStart, timeEnd) ? gravityfield->radialGradient(time, point) : 0;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double GravityfieldInInterval::field(const Time &time, const Vector3d &point, const Kernel &kernel2) const
{
  try
  {
    return time.isInInterval(timeStart, timeEnd) ? gravityfield->field(time, point, kernel2) : 0;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector3d GravityfieldInInterval::gravity(const Time &time, const Vector3d &point) const
{
  try
  {
    return time.isInInterval(timeStart, timeEnd) ? gravityfield->gravity(time, point) : Vector3d();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Tensor3d GravityfieldInInterval::gravityGradient(const Time &time, const Vector3d &point) const
{
  try
  {
    return time.isInInterval(timeStart, timeEnd) ? gravityfield->gravityGradient(time, point) : Tensor3d();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector3d GravityfieldInInterval::deformation(const Time &time, const Vector3d &point, Double gravity, const Vector &hn, const Vector &ln) const
{
  try
  {
    return time.isInInterval(timeStart, timeEnd) ? gravityfield->deformation(time, point, gravity, hn, ln) : Vector3d();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GravityfieldInInterval::deformation(const std::vector<Time> &time, const std::vector<Vector3d> &point, const std::vector<Double> &gravity,
                                         const Vector &hn, const Vector &ln, std::vector<std::vector<Vector3d>> &disp) const
{
  try
  {
    std::vector<Time> time2;
    std::copy_if(time.begin(), time.end(), std::back_inserter(time2), [&](const Time &x) {return ((x >= timeStart) && (x < timeEnd));});
    if(time2.size() == 0)
      return;

    std::vector<std::vector<Vector3d>> disp2(point.size(), std::vector<Vector3d>(time2.size()));
    gravityfield->deformation(time2, point, gravity, hn, ln, disp2);

    auto iter = time.begin();
    for(UInt k=0; k<time2.size(); k++)
    {
      const UInt index = std::distance(time.begin(), std::find(iter++, time.end(), time2.at(k)));
      for(UInt i=0; i<point.size(); i++)
        disp.at(i).at(index) = disp2.at(i).at(k);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

SphericalHarmonics GravityfieldInInterval::sphericalHarmonics(const Time &time, UInt maxDegree, UInt minDegree, Double GM, Double R) const
{
  try
  {
    return time.isInInterval(timeStart, timeEnd) ? gravityfield->sphericalHarmonics(time, maxDegree, minDegree, GM, R) : SphericalHarmonics().get(maxDegree, minDegree, GM, R);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GravityfieldInInterval::variance(const Time &time, const std::vector<Vector3d> &point, const Kernel &kernel, Matrix &D) const
{
  try
  {
    if(time.isInInterval(timeStart, timeEnd))
      D += gravityfield->variance(time, point, kernel);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
