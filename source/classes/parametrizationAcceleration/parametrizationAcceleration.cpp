/***********************************************/
/**
* @file parametrizationAcceleration.cpp
*
* @brief Orbit force parameters.
*
* @author Torsten Mayer-Guerr
* @date 2014-03-18
*
*/
/***********************************************/

#define DOCSTRING_ParametrizationAcceleration

#include "base/import.h"
#include "files/fileSatelliteModel.h"
#include "config/configRegister.h"
#include "classes/parametrizationAcceleration/parametrizationAccelerationPerRevolution.h"
#include "classes/parametrizationAcceleration/parametrizationAccelerationAccBias.h"
#include "classes/parametrizationAcceleration/parametrizationAccelerationAccScaleFactors.h"
#include "classes/parametrizationAcceleration/parametrizationAccelerationGnssSolarRadiation.h"
#include "classes/parametrizationAcceleration/parametrizationAccelerationThermosphericDensity.h"
#include "classes/parametrizationAcceleration/parametrizationAccelerationModelScale.h"
#include "classes/parametrizationAcceleration/parametrizationAcceleration.h"

/***********************************************/

GROOPS_REGISTER_CLASS(ParametrizationAcceleration, "parametrizationAccelerationType",
                      ParametrizationAccelerationPerRevolution,
                      ParametrizationAccelerationAccBias,
                      ParametrizationAccelerationAccScaleFactors,
                      ParametrizationAccelerationGnssSolarRadiation,
                      ParametrizationAccelerationThermosphericDensity,
                      ParametrizationAccelerationModelScale)

GROOPS_RENAMED_CLASS(parameterSatelliteType, parametrizationAccelerationType, date2time(2020, 6, 3))

GROOPS_READCONFIG_UNBOUNDED_CLASS(ParametrizationAcceleration, "parametrizationAccelerationType")

/***********************************************/

ParametrizationAcceleration::ParametrizationAcceleration(Config &config, const std::string &name)
{
  try
  {
    std::string type;
    while(readConfigChoice(config, name, type, Config::OPTIONAL, "", "orbit force parameters"))
    {
      if(readConfigChoiceElement(config, "perRevolution",             type, "oscillation per revoultion"))
        parameter.push_back(new ParametrizationAccelerationPerRevolution(config));
      if(readConfigChoiceElement(config, "accBias"  ,                 type, "accelerometer bias"))
        parameter.push_back(new ParametrizationAccelerationAccBias(config));
      if(readConfigChoiceElement(config, "accelerometerScaleFactors", type, "accelerometer scale factors"))
        parameter.push_back(new ParametrizationAccelerationAccScaleFactors(config));
      if(readConfigChoiceElement(config, "gnssSolarRadiation",        type, "GNSS solar radiation pressure model"))
        parameter.push_back(new ParametrizationAccelerationGnssSolarRadiation(config));
      if(readConfigChoiceElement(config, "thermosphericDensity",      type, "thermospheric density along the orbit"))
        parameter.push_back(new ParametrizationAccelerationThermosphericDensity(config));
      if(readConfigChoiceElement(config, "modelScale",                type, "model scale factor."))
        parameter.push_back(new ParametrizationAccelerationModelScale(config));
      endChoice(config);
      if(isCreateSchema(config))
        return;
    }

    computeIndicies();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

ParametrizationAcceleration::~ParametrizationAcceleration()
{
  for(UInt i=0; i<parameter.size(); i++)
    delete parameter.at(i);
}

/***********************************************/

void ParametrizationAcceleration::computeIndicies()
{
  try
  {
    parameterCountA = 0;
    parameterCountB = 0;
    indexA.resize(parameter.size());
    indexB.resize(parameter.size());
    for(UInt i=0; i<parameter.size(); i++)
    {
      indexA.at(i) = parameterCountA;
      indexB.at(i) = parameterCountB;
      if(parameter.at(i)->isPerArc())
        parameterCountB += parameter.at(i)->parameterCount();
      else
        parameterCountA += parameter.at(i)->parameterCount();
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool ParametrizationAcceleration::setInterval(const Time &timeStart, const Time &timeEnd)
{
  try
  {
    Bool change = FALSE;
    for(UInt i=0; i<parameter.size(); i++)
      if(!parameter.at(i)->isPerArc())
        change = parameter.at(i)->setInterval(timeStart, timeEnd) || change;
    if(change)
      computeIndicies();
    return change;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool ParametrizationAcceleration::setIntervalArc(const Time &timeStart, const Time &timeEnd)
{
  try
  {
    Bool change = FALSE;
    for(UInt i=0; i<parameter.size(); i++)
      if(parameter.at(i)->isPerArc())
        change = parameter.at(i)->setInterval(timeStart, timeEnd) || change;
    if(change)
      computeIndicies();
    return change;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ParametrizationAcceleration::parameterName(std::vector<ParameterName> &name) const
{
  for(UInt i=0; i<parameter.size(); i++)
    if(!parameter.at(i)->isPerArc())
      parameter.at(i)->parameterName(name);
}

/***********************************************/

void ParametrizationAcceleration::parameterNameArc(std::vector<ParameterName> &name) const
{
  for(UInt i=0; i<parameter.size(); i++)
    if(parameter.at(i)->isPerArc())
      parameter.at(i)->parameterName(name);
}

/***********************************************/

void ParametrizationAcceleration::compute(SatelliteModelPtr satellite, const Time &time, const Vector3d &position, const Vector3d &velocity,
                                 const Rotary3d &rotSat, const Rotary3d &rotEarth, EphemeridesPtr ephemerides, MatrixSliceRef A, MatrixSliceRef B)
{
  try
  {
    if(satellite)
    {
      Vector3d positionSun;
      if(ephemerides)
        positionSun = ephemerides->position(time, Ephemerides::SUN);
      satellite->changeState(time, position, velocity, positionSun, rotSat, rotEarth);
    }

    for(UInt i=0; i<parameter.size(); i++)
      parameter.at(i)->compute(satellite, time, position, velocity, rotSat, rotEarth, ephemerides,
                               (parameter.at(i)->isPerArc() ? B.slice(0, indexB.at(i), 3, parameter.at(i)->parameterCount()) :
                                                              A.slice(0, indexA.at(i), 3, parameter.at(i)->parameterCount())));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
