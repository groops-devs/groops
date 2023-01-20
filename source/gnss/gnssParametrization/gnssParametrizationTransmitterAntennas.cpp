/***********************************************/
/**
* @file gnssParametrizationTransmitterAntennas.cpp
*
* @brief Antenna center variations.
* @see GnssParametrization
*
* @author Torsten Mayer-Guerr
* @date 2021-01-23
*
*/
/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "classes/platformSelector/platformSelector.h"
#include "classes/parametrizationGnssAntenna/parametrizationGnssAntenna.h"
#include "gnss/gnss.h"
#include "gnss/gnssParametrization/gnssParametrizationTransmitterAntennas.h"

/***********************************************/

GnssParametrizationTransmitterAntennas::GnssParametrizationTransmitterAntennas(Config &config)
{
  try
  {
    readConfig(config, "name",                    name,                Config::OPTIONAL, "parameter.transmitterAntenna", "used for parameter selection");
    readConfig(config, "selectTransmitters",      selectTransmitters,  Config::MUSTSET,  "",  "");
    readConfig(config, "antennaCenterVariations", parametrization,     Config::MUSTSET,  "",  "estimate antenna center variations");
    readConfig(config, "patternTypes",            typesPattern,        Config::OPTIONAL, "",  "gnssType for each pattern (first match is used)");
    readConfig(config, "addNonMatchingTypes",     addNonMatchingTypes, Config::DEFAULT,  "1", "add patterns for additional observed gnssTypes that don't match any of the above");
    readConfig(config, "groupAntennas",           ignoreSerial,        Config::DEFAULT,  "0", "common ACVs for same antenna build types (ignores antenna serial number)");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationTransmitterAntennas::init(Gnss *gnss, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    this->gnss = gnss;
    transmitter2antenna.resize(gnss->transmitters.size(), NULLINDEX);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/


void GnssParametrizationTransmitterAntennas::initParameter(GnssNormalEquationInfo &normalEquationInfo)
{
  try
  {
    index.clear();
    types.clear();
    std::fill(transmitter2antenna.begin(), transmitter2antenna.end(), NULLINDEX);
    if(!isEnabled(normalEquationInfo, name) || normalEquationInfo.isEachReceiverSeparately)
      return;

    std::vector<std::string> antennaNames;
    auto selectedTransmitters = gnss->selectTransmitters(selectTransmitters);
    for(auto trans : gnss->transmitters)
      if(trans->useable() && selectedTransmitters.at(trans->idTrans()))
      {
        Bool found = FALSE;
        for(const auto &instrument : trans->platform.equipments)
        {
          auto antenna = std::dynamic_pointer_cast<PlatformGnssAntenna>(instrument);
          if(antenna && (antenna->timeEnd > gnss->times.front()) && (antenna->timeStart <= gnss->times.back()))
          {
            if(found)
            {
              logWarning<<trans->name()<<" has multiple antennas in interval -> only first one is used for parametrization/estimation"<<Log::endl;
              continue;
            }
            found = TRUE;

            // not already in list?
            const std::string name = GnssAntennaDefinition::str(antenna->name, ignoreSerial ? ""s : antenna->serial, antenna->radome);
            const UInt idAnt = std::distance(antennaNames.begin(), std::find(antennaNames.begin(), antennaNames.end(), name));
            if(idAnt >= antennaNames.size())
            {
              antennaNames.push_back(name);
              types.resize(antennaNames.size());
            }
            transmitter2antenna.at(trans->idTrans()) = idAnt;

            // add observed types
            for(const auto &typesRecv : gnss->typesRecvTrans)
              for(GnssType typeObs : GnssType::replaceCompositeSignals(typesRecv.at(trans->idTrans())))
                if(!typeObs.isInList(types.at(idAnt)))
                {
                  UInt idx;
                  if(typeObs.isInList(typesPattern, idx))
                    types.at(idAnt).push_back(typesPattern.at(idx));
                  else if(addNonMatchingTypes)
                    types.at(idAnt).push_back(typeObs);
                }
          } // for(antenna)
        }
      } // for(trans)

    // setup parameters
    // ----------------
    UInt countPara = 0;
    std::vector<ParameterName> baseNames;
    parametrization->parameterName(baseNames);
    index.resize(antennaNames.size());
    for(UInt idAnt=0; idAnt<antennaNames.size(); idAnt++)
      if(types.at(idAnt).size())
      {
        std::sort(types.at(idAnt).begin(), types.at(idAnt).end());
        for(GnssType type : types.at(idAnt))
        {
          std::string typeStr = "." + type.str();
          std::vector<ParameterName> parameterNames;
          for(const auto &base : baseNames)
            parameterNames.push_back(ParameterName(antennaNames.at(idAnt), base.type + typeStr, base.temporal, base.interval));
          index.at(idAnt).push_back(normalEquationInfo.parameterNamesOther(parameterNames));
          countPara += parameterNames.size();
        }
      }
    if(countPara)
      logInfo<<countPara%"%9i transmitter antenna parameters"s<<Log::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationTransmitterAntennas::designMatrix(const GnssNormalEquationInfo &/*normalEquationInfo*/, const GnssObservationEquation &eqn, GnssDesignMatrix &A) const
{
  try
  {
    const UInt idAnt = transmitter2antenna.at(eqn.transmitter->idTrans());
    if((idAnt >= index.size()) || !index.at(idAnt).size())
      return;

    Matrix B(eqn.l.rows(), types.at(idAnt).size());
    std::vector<Bool> used(types.at(idAnt).size(), FALSE);
    UInt idPattern;
    for(UInt idType=0; idType<eqn.typesTransmitted.size(); idType++)
      if(eqn.typesTransmitted.at(idType).isInList(types.at(idAnt), idPattern))
      {
        axpy(1., eqn.A.column(GnssObservationEquation::idxUnit + eqn.types.size() + idType), B.column(idPattern));
        used.at(idPattern) = TRUE;
      }
    if(!std::any_of(used.begin(), used.end(), [](Bool x){return x;}))
      return;

    const Matrix Acv = parametrization->designMatrix(eqn.azimutTrans, eqn.elevationTrans);
    for(UInt idPattern=0; idPattern<types.at(idAnt).size(); idPattern++)
      if(used.at(idPattern))
        matMult(1., B.column(idPattern), Acv, A.column(index.at(idAnt).at(idPattern)));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
