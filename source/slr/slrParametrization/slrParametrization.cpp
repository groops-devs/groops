/***********************************************/
/**
* @file slrParametrization.cpp
*
* @brief Parametrization of SLR observations.
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#define DOCSTRING_SlrParametrization

#include "base/import.h"
#include "config/configRegister.h"
#include "slr/slrParametrization/slrParametrizationTroposphere.h"
#include "slr/slrParametrization/slrParametrizationDynamicOrbits.h"
#include "slr/slrParametrization/slrParametrizationGravityField.h"
#include "slr/slrParametrization/slrParametrizationStaticPositions.h"
#include "slr/slrParametrization/slrParametrizationEarthRotation.h"
#include "slr/slrParametrization/slrParametrizationRangeBiasSatelliteApriori.h"
#include "slr/slrParametrization/slrParametrizationRangeBiasSatellite.h"
#include "slr/slrParametrization/slrParametrizationRangeBiasStationApriori.h"
#include "slr/slrParametrization/slrParametrizationRangeBiasStation.h"
#include "slr/slrParametrization/slrParametrizationRangeBiasStationSatellite.h"
#include "slr/slrParametrization/slrParametrizationRangeBiasStationSatelliteApriori.h"
#include "slr/slrParametrization/slrParametrizationTimeBiasApriori.h"
#include "slr/slrParametrization/slrParametrizationTimeBias.h"
#include "slr/slrParametrization/slrParametrizationConstraints.h"
#include "slr/slrParametrization/slrParametrizationGroup.h"
#include "slr/slrParametrization/slrParametrization.h"

/***********************************************/

GROOPS_REGISTER_CLASS(SlrParametrization, "slrParametrizationType",
                      SlrParametrizationTroposphere,
                      SlrParametrizationDynamicOrbits,
                      SlrParametrizationGravityField,
                      SlrParametrizationStaticPositions,
                      SlrParametrizationEarthRotation,
                      SlrParametrizationRangeBiasStationApriori,
                      SlrParametrizationRangeBiasStation,
                      SlrParametrizationRangeBiasSatelliteApriori,
                      SlrParametrizationRangeBiasSatellite,
                      SlrParametrizationRangeBiasStationSatelliteApriori,
                      SlrParametrizationRangeBiasStationSatellite,
                      SlrParametrizationTimeBiasApriori,
                      SlrParametrizationTimeBias,
                      SlrParametrizationConstraints,
                      SlrParametrizationGroup)

GROOPS_READCONFIG_UNBOUNDED_CLASS(SlrParametrization, "slrParametrizationType")

/***********************************************/

SlrParametrization::SlrParametrization(Config &config, const std::string &name)
{
  try
  {
    std::string type;
    while(readConfigChoice(config, name, type, Config::OPTIONAL, "", "parametrization of SLR observations"))
    {
      if(readConfigChoiceElement(config, "troposphere",                      type, "tropospheric delays"))
        base.push_back(new SlrParametrizationTroposphere(config));
      if(readConfigChoiceElement(config, "dynamicOrbits",                    type, "satellite orbits by variational equations"))
        base.push_back(new SlrParametrizationDynamicOrbits(config));
      if(readConfigChoiceElement(config, "gravityField",                     type, "gravity field of the Earth"))
        base.push_back(new SlrParametrizationGravityField(config));
      if(readConfigChoiceElement(config, "staticPositions",                  type, "static positions with no-net constraints"))
        base.push_back(new SlrParametrizationStaticPositions(config));
      if(readConfigChoiceElement(config, "earthRotation",                    type, "Earth rotation"))
        base.push_back(new SlrParametrizationEarthRotation(config));
      if(readConfigChoiceElement(config, "rangeBiasStationApriori",          type, "apriori range bias from file"))
        base.push_back(new SlrParametrizationRangeBiasStationApriori(config));
      if(readConfigChoiceElement(config, "rangeBiasStation",                 type, "range bias"))
        base.push_back(new SlrParametrizationRangeBiasStation(config));
      if(readConfigChoiceElement(config, "rangeBiasSatelliteApriori",        type, "apriori range bias from file"))
        base.push_back(new SlrParametrizationRangeBiasSatelliteApriori(config));
      if(readConfigChoiceElement(config, "rangeBiasSatellite",               type, "range bias"))
        base.push_back(new SlrParametrizationRangeBiasSatellite(config));
      if(readConfigChoiceElement(config, "rangeBiasStationSatelliteApriori", type, "apriori range bias from file"))
        base.push_back(new SlrParametrizationRangeBiasStationSatelliteApriori(config));
      if(readConfigChoiceElement(config, "rangeBiasStationSatellite",        type, "range bias"))
        base.push_back(new SlrParametrizationRangeBiasStationSatellite(config));
      if(readConfigChoiceElement(config, "timeBiasApriori",                  type, "apriori time bias at station from file"))
        base.push_back(new SlrParametrizationTimeBiasApriori(config));
      if(readConfigChoiceElement(config, "timeBias",                         type, "time bias at station"))
        base.push_back(new SlrParametrizationTimeBias(config));
      if(readConfigChoiceElement(config, "constraints",                      type, "parameter constraints"))
        base.push_back(new SlrParametrizationConstraints(config));
      if(readConfigChoiceElement(config, "group",                            type, "grouping parametrizations"))
        base.push_back(new SlrParametrizationGroup(config));
      endChoice(config);
      if(isCreateSchema(config))
        return;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

SlrParametrization::~SlrParametrization()
{
  for(auto b : base)
    delete b;
}

/***********************************************/

std::vector<const SlrParametrizationGravityField*> SlrParametrization::getParametrizationGravity() const
{
  try
  {
    std::vector<const SlrParametrizationGravityField*> params;
    for(auto b : base)
      b->getParametrizationGravity(params);
    return params;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SlrParametrization::init(Slr *slr, const std::vector<const SlrParametrizationGravityField*> &paramGravityField)
{
  try
  {
    for(auto b : base)
      b->init(slr, paramGravityField);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SlrParametrization::observationCorrections(SlrObservationEquation &eqn) const
{
  try
  {
    for(auto b : base)
      b->observationCorrections(eqn);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SlrParametrization::initParameter(SlrNormalEquationInfo &normalEquationInfo)
{
  try
  {
    for(auto b : base)
      b->initParameter(normalEquationInfo);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector SlrParametrization::aprioriParameter(const SlrNormalEquationInfo &normalEquationInfo) const
{
  try
  {
    Vector x0(normalEquationInfo.parameterCount());
    for(auto b : base)
      b->aprioriParameter(normalEquationInfo, x0);
    return x0;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SlrParametrization::designMatrix(const SlrNormalEquationInfo &normalEquationInfo, const SlrObservationEquation &eqn, SlrDesignMatrix &A) const
{
  try
  {
    for(auto b : base)
      b->designMatrix(normalEquationInfo, eqn, A);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SlrParametrization::constraints(const SlrNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const
{
  try
  {
    for(auto b : base)
      b->constraints(normalEquationInfo, normals, n, lPl, obsCount);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double SlrParametrization::updateParameter(const SlrNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz)
{
  try
  {
    Double change = 0;
    for(auto b : base)
      change = std::max(change, b->updateParameter(normalEquationInfo, x, Wz));
    return change;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SlrParametrization::updateCovariance(const SlrNormalEquationInfo &normalEquationInfo, const MatrixDistributed &covariance)
{
  try
  {
    for(auto b : base)
      b->updateCovariance(normalEquationInfo, covariance);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SlrParametrization::writeResults(const SlrNormalEquationInfo &normalEquationInfo, const std::string &suffix) const
{
  try
  {
    for(auto b : base)
      b->writeResults(normalEquationInfo, suffix);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

Bool SlrParametrizationBase::isEnabled(const SlrNormalEquationInfo &normalEquationInfo, const std::string &name)
{
  try
  {
    Bool enabled = TRUE;
    for(const auto &pattern : normalEquationInfo.enableParametrizations)
      if(std::regex_match(name, pattern.first))
        enabled = pattern.second;
    return enabled;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
