/***********************************************/
/**
* @file parametrizationSatelliteTracking.cpp
*
* @brief SST parameter.
*
* @author Beate Klinger
* @date 2015-04-27
*
*/
/***********************************************/

#define DOCSTRING_ParametrizationSatelliteTracking

#include "base/import.h"
#include "base/parameterName.h"
#include "base/sphericalHarmonics.h"
#include "config/configRegister.h"
#include "classes/parametrizationSatelliteTracking/parametrizationSatelliteTrackingAntennaCenter.h"
#include "classes/parametrizationSatelliteTracking/parametrizationSatelliteTrackingBias.h"
#include "classes/parametrizationSatelliteTracking/parametrizationSatelliteTrackingScale.h"
#include "classes/parametrizationSatelliteTracking/parametrizationSatelliteTrackingTimeBias.h"
#include "classes/parametrizationSatelliteTracking/parametrizationSatelliteTrackingScaleModel.h"
#include "classes/parametrizationSatelliteTracking/parametrizationSatelliteTrackingSpecialEffect.h"
#include "classes/parametrizationSatelliteTracking/parametrizationSatelliteTracking.h"

/***********************************************/

GROOPS_REGISTER_CLASS(ParametrizationSatelliteTracking, "parametrizationSatelliteTrackingType",
                      ParametrizationSatelliteTrackingAntennaCenter,
                      ParametrizationSatelliteTrackingBias,
                      ParametrizationSatelliteTrackingScale,
                      ParametrizationSatelliteTrackingTimeBias,
                      ParametrizationSatelliteTrackingScaleModel,
                      ParametrizationSatelliteTrackingSpecialEffect)

GROOPS_RENAMED_CLASS(parameterSatelliteTrackingType, parametrizationSatelliteTrackingType, date2time(2020, 6, 3))

GROOPS_READCONFIG_UNBOUNDED_CLASS(ParametrizationSatelliteTracking, "parametrizationSatelliteTrackingType")

/***********************************************/

ParametrizationSatelliteTracking::ParametrizationSatelliteTracking(Config &config, const std::string &name)
{
  try
  {
    std::string type;
    while(readConfigChoice(config, name, type, Config::OPTIONAL, "", "SST parameters"))
    {
      if(readConfigChoiceElement(config, "antennaCenter",            type, "antenna center"))
        parameter.push_back(new ParametrizationSatelliteTrackingAntennaCenter(config));
      if(readConfigChoiceElement(config, "bias",                     type, "SST bias"))
        parameter.push_back(new ParametrizationSatelliteTrackingBias(config));
      if(readConfigChoiceElement(config, "scale",                    type, "SST scale"))
        parameter.push_back(new ParametrizationSatelliteTrackingScale(config));
      if(readConfigChoiceElement(config, "timeBias",                 type, "SST time bias"))
        parameter.push_back(new ParametrizationSatelliteTrackingTimeBias(config));
      if(readConfigChoiceElement(config, "scaleModel",               type, "scale factors for model from file"))
        parameter.push_back(new ParametrizationSatelliteTrackingScaleModel(config));
      if(readConfigChoiceElement(config, "specialEffect",            type, "special effects from file"))
        parameter.push_back(new ParametrizationSatelliteTrackingSpecialEffect(config));
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

ParametrizationSatelliteTracking::~ParametrizationSatelliteTracking()
{
  for(UInt i=0; i<parameter.size(); i++)
    delete parameter.at(i);
}

/***********************************************/

void ParametrizationSatelliteTracking::computeIndicies()
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

void ParametrizationSatelliteTracking::setInterval(const Time &timeStart, const Time &timeEnd)
{
  try
  {
    for(UInt i=0; i<parameter.size(); i++)
      if(!parameter.at(i)->isPerArc())
        parameter.at(i)->setInterval(timeStart, timeEnd);
    computeIndicies();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ParametrizationSatelliteTracking::setIntervalArc(const Time &timeStart, const Time &timeEnd)
{
  try
  {
    for(UInt i=0; i<parameter.size(); i++)
      if(parameter.at(i)->isPerArc())
        parameter.at(i)->setInterval(timeStart, timeEnd);
    computeIndicies();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ParametrizationSatelliteTracking::parameterName(std::vector<ParameterName> &name) const
{
  for(UInt i=0; i<parameter.size(); i++)
    if(!parameter.at(i)->isPerArc())
      parameter.at(i)->parameterName(name);
}

/***********************************************/

void ParametrizationSatelliteTracking::parameterNameArc(std::vector<ParameterName> &name) const
{
  for(UInt i=0; i<parameter.size(); i++)
    if(parameter.at(i)->isPerArc())
      parameter.at(i)->parameterName(name);
}

/***********************************************/

void ParametrizationSatelliteTracking::compute(UInt sstType, const std::vector<Time> &time, const Vector &sst0,
                                         const Vector &position1, const Vector &position2, const Vector &velocity1, const Vector &velocity2,
                                         const std::vector<Rotary3d> &rotSat1, const std::vector<Rotary3d> &rotSat2, MatrixSliceRef A, MatrixSliceRef B)
{
  try
  {
    for(UInt i=0; i<parameter.size(); i++)
      parameter.at(i)->compute(sstType, time, sst0, position1, position2, velocity1, velocity2, rotSat1, rotSat2,
                               (parameter.at(i)->isPerArc() ? B.slice(0, indexB.at(i), time.size(), parameter.at(i)->parameterCount()) :
                                                              A.slice(0, indexA.at(i), time.size(), parameter.at(i)->parameterCount())));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
