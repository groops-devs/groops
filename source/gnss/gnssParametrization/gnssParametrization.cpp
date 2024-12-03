/***********************************************/
/**
* @file gnssParametrization.cpp
*
* @brief Parametrization of GNSS observations.
*
* @author Torsten Mayer-Guerr
* @date 2021-01-23
*
*/
/***********************************************/

#define DOCSTRING_GnssParametrization

#include "base/import.h"
#include "config/configRegister.h"
#include "gnss/gnssParametrization/gnssParametrizationIonosphereSTEC.h"
#include "gnss/gnssParametrization/gnssParametrizationIonosphereVTEC.h"
#include "gnss/gnssParametrization/gnssParametrizationIonosphereMap.h"
#include "gnss/gnssParametrization/gnssParametrizationClocks.h"
#include "gnss/gnssParametrization/gnssParametrizationClocksModel.h"
#include "gnss/gnssParametrization/gnssParametrizationSignalBiases.h"
#include "gnss/gnssParametrization/gnssParametrizationAmbiguities.h"
#include "gnss/gnssParametrization/gnssParametrizationCodeBiases.h"
#include "gnss/gnssParametrization/gnssParametrizationTecBiases.h"
#include "gnss/gnssParametrization/gnssParametrizationTemporalBias.h"
#include "gnss/gnssParametrization/gnssParametrizationStaticPositions.h"
#include "gnss/gnssParametrization/gnssParametrizationKinematicPositions.h"
#include "gnss/gnssParametrization/gnssParametrizationLeoDynamicOrbits.h"
#include "gnss/gnssParametrization/gnssParametrizationTransmitterDynamicOrbits.h"
#include "gnss/gnssParametrization/gnssParametrizationTroposphere.h"
#include "gnss/gnssParametrization/gnssParametrizationEarthRotation.h"
#include "gnss/gnssParametrization/gnssParametrizationReceiverAntennas.h"
#include "gnss/gnssParametrization/gnssParametrizationTransmitterAntennas.h"
#include "gnss/gnssParametrization/gnssParametrizationConstraints.h"
#include "gnss/gnssParametrization/gnssParametrizationGroup.h"
#include "gnss/gnssParametrization/gnssParametrization.h"

/***********************************************/

GROOPS_REGISTER_CLASS(GnssParametrization, "gnssParametrizationType",
                      GnssParametrizationIonosphereSTEC,
                      GnssParametrizationIonosphereVTEC,
                      GnssParametrizationIonosphereMap,
                      GnssParametrizationClocks,
                      GnssParametrizationClocksModel,
                      GnssParametrizationSignalBiases,
                      GnssParametrizationAmbiguities,
                      GnssParametrizationCodeBiases,
                      GnssParametrizationTecBiases,
                      GnssParametrizationTemporalBias,
                      GnssParametrizationStaticPositions,
                      GnssParametrizationKinematicPositions,
                      GnssParametrizationLeoDynamicOrbits,
                      GnssParametrizationTransmitterDynamicOrbits,
                      GnssParametrizationTroposphere,
                      GnssParametrizationEarthRotation,
                      GnssParametrizationReceiverAntennas,
                      GnssParametrizationTransmitterAntennas,
                      GnssParametrizationConstraints,
                      GnssParametrizationGroup)

GROOPS_READCONFIG_UNBOUNDED_CLASS(GnssParametrization, "gnssParametrizationType")

/***********************************************/

GnssParametrization::GnssParametrization(Config &config, const std::string &name)
{
  try
  {
    std::string type;
    while(readConfigChoice(config, name, type, Config::OPTIONAL, "", "parametrization of GNSS observations"))
    {
      if(readConfigChoiceElement(config, "ionosphereSTEC",           type, "ionospheric slant delays"))
        base.push_back(new GnssParametrizationIonosphereSTEC(config));
      if(readConfigChoiceElement(config, "ionosphereVTEC",           type, "ionospheric vertical delays"))
        base.push_back(new GnssParametrizationIonosphereVTEC(config));
      if(readConfigChoiceElement(config, "ionosphereMap",            type, "ionospheric vertical delays"))
        base.push_back(new GnssParametrizationIonosphereMap(config));
      if(readConfigChoiceElement(config, "clocks",                   type, "clock errors"))
        base.push_back(new GnssParametrizationClocks(config));
      if(readConfigChoiceElement(config, "clocksModel",              type, "clock errors with AR1 process model"))
        base.push_back(new GnssParametrizationClocksModel(config));
      if(readConfigChoiceElement(config, "signalBiases",             type, "apriori values of signal biases (code/phase)"))
        base.push_back(new GnssParametrizationSignalBiases(config));
      if(readConfigChoiceElement(config, "ambiguities",              type, "integer and float ambiguities"))
        base.push_back(new GnssParametrizationAmbiguities(config));
      if(readConfigChoiceElement(config, "codeBiases",               type, "biases of code signals"))
        base.push_back(new GnssParametrizationCodeBiases(config));
      if(readConfigChoiceElement(config, "tecBiases",                type, "biases of signals which changes the estimated TEC"))
        base.push_back(new GnssParametrizationTecBiases(config));
      if(readConfigChoiceElement(config, "temporalBias",             type, "temporal changing signal bias"))
        base.push_back(new GnssParametrizationTemporalBias(config));
      if(readConfigChoiceElement(config, "staticPositions",          type, "static positions with no-net constraints"))
        base.push_back(new GnssParametrizationStaticPositions(config));
      if(readConfigChoiceElement(config, "kinematicPositions",       type, "position each epoch"))
        base.push_back(new GnssParametrizationKinematicPositions(config));
      if(readConfigChoiceElement(config, "leoDynamicOrbits",          type, "LEO orbits by variational equations"))
        base.push_back(new GnssParametrizationLeoDynamicOrbits(config));
      if(readConfigChoiceElement(config, "transmitterDynamicOrbits", type, "GNSS satellite orbits by variational equations"))
        base.push_back(new GnssParametrizationTransmitterDynamicOrbits(config));
      if(readConfigChoiceElement(config, "troposphere",              type, "tropospheric delays"))
        base.push_back(new GnssParametrizationTroposphere(config));
      if(readConfigChoiceElement(config, "earthRotation",            type, "Earth rotation"))
        base.push_back(new GnssParametrizationEarthRotation(config));
      if(readConfigChoiceElement(config, "receiverAntennas",         type, "antenna center variations"))
        base.push_back(new GnssParametrizationReceiverAntennas(config));
      if(readConfigChoiceElement(config, "transmitterAntennas",      type, "antenna center variations"))
        base.push_back(new GnssParametrizationTransmitterAntennas(config));
      if(readConfigChoiceElement(config, "constraints",              type, "parameter constraints"))
        base.push_back(new GnssParametrizationConstraints(config));
      if(readConfigChoiceElement(config, "group",                    type, "grouping parametrizations"))
        base.push_back(new GnssParametrizationGroup(config));
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

GnssParametrization::~GnssParametrization()
{
  for(auto b : base)
    delete b;
}

/***********************************************/

void GnssParametrization::init(Gnss *gnss, Parallel::CommunicatorPtr comm)
{
  try
  {
    for(auto b : base)
      b->init(gnss, comm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrization::observationCorrections(GnssObservationEquation &eqn) const
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

void GnssParametrization::requirements(GnssNormalEquationInfo &normalEquationInfo,
                                       std::vector<UInt> &transCount, std::vector<UInt> &transCountEpoch,
                                       std::vector<UInt> &recvCount, std::vector<UInt> &recvCountEpoch)
{
  try
  {
    for(auto b : base)
      b->requirements(normalEquationInfo, transCount, transCountEpoch, recvCount, recvCountEpoch);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrization::initParameter(GnssNormalEquationInfo &normalEquationInfo)
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

Vector GnssParametrization::aprioriParameter(const GnssNormalEquationInfo &normalEquationInfo) const
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

void GnssParametrization::designMatrix(const GnssNormalEquationInfo &normalEquationInfo, const GnssObservationEquation &eqn, GnssDesignMatrix &A) const
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

void GnssParametrization::constraintsEpoch(const GnssNormalEquationInfo &normalEquationInfo, UInt idEpoch, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const
{
  try
  {
    for(auto b : base)
      b->constraintsEpoch(normalEquationInfo, idEpoch, normals, n, lPl, obsCount);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrization::constraints(const GnssNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const
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

Double GnssParametrization::ambiguityResolve(const GnssNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount,
                                             const std::vector<Byte> &selectedTransmitters, const std::vector<Byte> &selectedReceivers,
                                             const std::function<Vector(const_MatrixSliceRef xFloat, MatrixSliceRef W, const_MatrixSliceRef d, Vector &xInt, Double &sigma)> &searchInteger)
{
  try
  {
    Double sigmaFloat = 0;
    for(auto b : base)
      sigmaFloat = std::max(sigmaFloat, b->ambiguityResolve(normalEquationInfo, normals, n, lPl, obsCount,
                                                            selectedTransmitters, selectedReceivers, searchInteger));
    return sigmaFloat;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}


/***********************************************/

Double GnssParametrization::updateParameter(const GnssNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz)
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

void GnssParametrization::updateCovariance(const GnssNormalEquationInfo &normalEquationInfo, const MatrixDistributed &covariance)
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

void GnssParametrization::writeResults(const GnssNormalEquationInfo &normalEquationInfo, const std::string &suffix) const
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

Bool GnssParametrizationBase::isEnabled(const GnssNormalEquationInfo &normalEquationInfo, const std::string &name)
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
