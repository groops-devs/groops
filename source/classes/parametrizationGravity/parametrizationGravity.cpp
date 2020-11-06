/***********************************************/
/**
* @file parametrizationGravity.cpp
*
* @brief Parametrization of the gravity field.
*
* @author Torsten Mayer-Guerr
* @date 2001-08-25
*
*/
/***********************************************/

#define DOCSTRING_ParametrizationGravity

#include "base/import.h"
#include "base/sphericalHarmonics.h"
#include "config/configRegister.h"
#include "classes/kernel/kernel.h"
#include "classes/parametrizationGravity/parametrizationGravitySphericalHarmonics.h"
#include "classes/parametrizationGravity/parametrizationGravityRadialBasis.h"
#include "classes/parametrizationGravity/parametrizationGravityTemporal.h"
#include "classes/parametrizationGravity/parametrizationGravityEarthquakeOscillation.h"
#include "classes/parametrizationGravity/parametrizationGravity.h"

/***********************************************/

GROOPS_REGISTER_CLASS(ParametrizationGravity, "parametrizationGravityType",
                      ParametrizationGravitySphericalHarmonics,
                      ParametrizationGravityRadialBasis,
                      ParametrizationGravityTemporal,
                      ParametrizationGravityEarthquakeOscillation)

GROOPS_RENAMED_CLASS(representationType, parametrizationGravityType, date2time(2020, 6, 3))

GROOPS_READCONFIG_UNBOUNDED_CLASS(ParametrizationGravity, "parametrizationGravityType")

/***********************************************/

ParametrizationGravity::ParametrizationGravity(Config &config, const std::string &name)
{
  try
  {
    std::string type;
    while(readConfigChoice(config, name, type, Config::OPTIONAL, "", "parametrization of the gravity field"))
    {
      if(readConfigChoiceElement(config, "sphericalHarmonics", type, "potential coefficents!"))
        parametrizations.push_back(new ParametrizationGravitySphericalHarmonics(config));
      if(readConfigChoiceElement(config, "radialBasis",        type, "harmonic radial basis functions"))
        parametrizations.push_back(new ParametrizationGravityRadialBasis(config));
      if(readConfigChoiceElement(config, "temporal",           type, "time variable gravity field"))
        parametrizations.push_back(new ParametrizationGravityTemporal(config));
      if(readConfigChoiceElement(config, "earthquakeOscillation", type, "earthquake oscillation parameters"))
        parametrizations.push_back(new ParametrizationGravityEarthquakeOscillation(config));
      endChoice(config);
      if(isCreateSchema(config))
        return;
    }

    computeIndices();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

ParametrizationGravity::~ParametrizationGravity()
{
  for(UInt i=0; i<parametrizations.size(); i++)
    delete parametrizations.at(i);
}

/***********************************************/
/***********************************************/

void ParametrizationGravity::computeIndices()
{
  try
  {
    _parameterCount = 0;
    index.resize(parametrizations.size());
    for(UInt i=0; i<parametrizations.size(); i++)
    {
      index.at(i) = _parameterCount;
      _parameterCount += parametrizations.at(i)->parameterCount();
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ParametrizationGravity::setInterval(const Time &timeStart, const Time &timeEnd)
{
  for(UInt i=0; i<parametrizations.size(); i++)
    parametrizations.at(i)->setInterval(timeStart, timeEnd);
  computeIndices();
}

/***********************************************/

void ParametrizationGravity::parameterName(std::vector<ParameterName> &name) const
{
  for(UInt i=0; i<parametrizations.size(); i++)
    parametrizations.at(i)->parameterName(name);
}

/***********************************************/

void ParametrizationGravity::field(const Time &time, const Vector3d &point, const Kernel &kernel, MatrixSliceRef A) const
{
  for(UInt i=0; i<parametrizations.size(); i++)
    parametrizations.at(i)->field(time, point, kernel, A.slice(0,index.at(i),1,parametrizations.at(i)->parameterCount()));
}

/***********************************************/

void ParametrizationGravity::potential(const Time &time, const Vector3d &point, MatrixSliceRef A) const
{
  for(UInt i=0; i<parametrizations.size(); i++)
    parametrizations.at(i)->potential(time, point, A.slice(0,index.at(i),1,parametrizations.at(i)->parameterCount()));
}

/***********************************************/

void ParametrizationGravity::radialGradient(const Time &time, const Vector3d &point, MatrixSliceRef A) const
{
  for(UInt i=0; i<parametrizations.size(); i++)
    parametrizations.at(i)->radialGradient(time, point, A.slice(0,index.at(i),1,parametrizations.at(i)->parameterCount()));
}

/***********************************************/

void ParametrizationGravity::gravity(const Time &time, const Vector3d &point, MatrixSliceRef A) const
{
  for(UInt i=0; i<parametrizations.size(); i++)
    parametrizations.at(i)->gravity(time, point, A.slice(0,index.at(i),3,parametrizations.at(i)->parameterCount()));
}

/***********************************************/

void ParametrizationGravity::gravityGradient(const Time &time, const Vector3d &point, MatrixSliceRef A) const
{
  for(UInt i=0; i<parametrizations.size(); i++)
    parametrizations.at(i)->gravityGradient(time, point, A.slice(0,index.at(i),6,parametrizations.at(i)->parameterCount()));
}

/***********************************************/

void ParametrizationGravity::deformation(const Time &time, const Vector3d &point, Double gravity, const Vector &hn, const Vector &ln, MatrixSliceRef A) const
{
  for(UInt i=0; i<parametrizations.size(); i++)
    parametrizations.at(i)->deformation(time, point, gravity, hn, ln, A.slice(0,index.at(i),3,parametrizations.at(i)->parameterCount()));
}

/***********************************************/

SphericalHarmonics ParametrizationGravity::sphericalHarmonics(const Time &time, const Vector &x, UInt maxDegree) const
{
  try
  {
    SphericalHarmonics harmonics;
    for(UInt i=0; i<parametrizations.size(); i++)
      harmonics += parametrizations.at(i)->sphericalHarmonics(time, x.slice(index.at(i),parametrizations.at(i)->parameterCount()), maxDegree);
    return harmonics;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

SphericalHarmonics ParametrizationGravity::sphericalHarmonics(const Time &time, const Vector &x, const Vector &sigma2x, UInt maxDegree) const
{
  try
  {
    SphericalHarmonics harmonics;
    for(UInt i=0; i<parametrizations.size(); i++)
      harmonics += parametrizations.at(i)->sphericalHarmonics(time, x.slice(index.at(i),parametrizations.at(i)->parameterCount()), sigma2x.slice(index.at(i),parametrizations.at(i)->parameterCount()), maxDegree);
    return harmonics;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
