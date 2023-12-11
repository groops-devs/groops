/***********************************************/
/**
 * @file gravityfieldGroup.cpp
 *
 * @brief Group.
 * @see Gravityfield
 *
 * @author Torsten Mayer-Guerr
 * @date 2023-11-13
 *
 */
/***********************************************/

#include "base/import.h"
#include "base/sphericalHarmonics.h"
#include "config/config.h"
#include "classes/kernel/kernel.h"
#include "classes/gravityfield/gravityfield.h"
#include "classes/gravityfield/gravityfieldGroup.h"

/***********************************************/

GravityfieldGroup::GravityfieldGroup(Config &config)
{
  try
  {
    readConfig(config, "gravityfield", gravityfield, Config::DEFAULT, "", "");
    readConfig(config, "factor",       factor,       Config::DEFAULT, "1.0", "the result is multiplied by this factor, set -1 to subtract the field");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double GravityfieldGroup::potential(const Time &time, const Vector3d &point) const
{
  return factor * gravityfield->potential(time, point);
}

/***********************************************/

Double GravityfieldGroup::radialGradient(const Time &time, const Vector3d &point) const
{
  return factor * gravityfield->radialGradient(time, point);
}

/***********************************************/

Double GravityfieldGroup::field(const Time &time, const Vector3d &point, const Kernel &kernel) const
{
  return factor * gravityfield->field(time, point, kernel);
}

/***********************************************/

Vector3d GravityfieldGroup::gravity(const Time &time, const Vector3d &point) const
{
  return gravityfield->gravity(time, point);
}

/***********************************************/

Tensor3d GravityfieldGroup::gravityGradient(const Time &time, const Vector3d &point) const
{
  return factor * gravityfield->gravityGradient(time, point);
}

/***********************************************/

Vector3d GravityfieldGroup::deformation(const Time &time, const Vector3d &point, Double gravity, const Vector &hn, const Vector &ln) const
{
  return factor * gravityfield->deformation(time, point, gravity, hn, ln);
}

/***********************************************/

void GravityfieldGroup::deformation(const std::vector<Time> &time, const std::vector<Vector3d> &point, const std::vector<Double> &gravity,
                                    const Vector &hn, const Vector &ln, std::vector<std::vector<Vector3d>> &disp) const
{
  gravityfield->deformation(time, point, gravity, hn, ln, disp);
  if(factor != 1.)
    for(auto &ds : disp)
      for(auto &d : ds)
        d *= factor;
}

/***********************************************/

void GravityfieldGroup::variance(const Time &time, const std::vector<Vector3d> &point, const Kernel &kernel, Matrix &D) const
{
  axpy(factor*factor, gravityfield->variance(time, point, kernel), D);
}

/***********************************************/

SphericalHarmonics GravityfieldGroup::sphericalHarmonics(const Time &time, UInt maxDegree, UInt minDegree, Double GM, Double R) const
{
  return factor * gravityfield->sphericalHarmonics(time, maxDegree, minDegree, GM, R);
}

/***********************************************/
