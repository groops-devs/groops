/***********************************************/
/**
* @file gnssParametrizationReceiverAntennas.cpp
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
#include "gnss/gnss.h"
#include "gnss/gnssTransceiverSelector/gnssTransceiverSelector.h"
#include "classes/parametrizationGnssAntenna/parametrizationGnssAntenna.h"
#include "gnss/gnssParametrization/gnssParametrizationReceiverAntennas.h"

/***********************************************/

GnssParametrizationReceiverAntennas::GnssParametrizationReceiverAntennas(Config &config)
{
  try
  {
    readConfig(config, "name",                    name,                Config::OPTIONAL, "parameter.receiverAntenna", "used for parameter selection");
    readConfig(config, "selectReceivers",         selectReceivers,     Config::MUSTSET,  "",  "");
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

void GnssParametrizationReceiverAntennas::init(Gnss *gnss, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    this->gnss = gnss;
    receiver2antenna.resize(gnss->receivers.size(), NULLINDEX);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/


void GnssParametrizationReceiverAntennas::initParameter(GnssNormalEquationInfo &normalEquationInfo)
{
  try
  {
    index.clear();
    types.clear();
    std::fill(receiver2antenna.begin(), receiver2antenna.end(), NULLINDEX);
    if(!isEnabled(normalEquationInfo, name))
      return;

    std::vector<std::string> antennaNames;
    auto selectedReceivers = selectReceivers->select(gnss->receivers);
    for(auto recv : gnss->receivers)
      if(recv->useable() && selectedReceivers.at(recv->idRecv()) && normalEquationInfo.estimateReceiver.at(recv->idRecv()))
      {
        Bool found = FALSE;
        for(const auto &antenna : recv->info.antenna)
          if((antenna.timeEnd > gnss->times.front()) && (antenna.timeStart <= gnss->times.back()))
          {
            if(found)
            {
              logWarning<<recv->name()<<" has multiple antennas in interval -> only first one is used for parametrization/estimation"<<Log::endl;
              continue;
            }
            found = TRUE;

            // not already in list?
            const std::string name = GnssAntennaDefinition::str(antenna.name, ignoreSerial ? ""s : antenna.serial, antenna.radome);
            const UInt idAnt = std::distance(antennaNames.begin(), std::find(antennaNames.begin(), antennaNames.end(), name));
            if(idAnt >= antennaNames.size())
            {
              antennaNames.push_back(name);
              types.resize(antennaNames.size());
            }
            receiver2antenna.at(recv->idRecv()) = idAnt;

            // add observed types
            for(const auto &typesRecv : gnss->typesRecvTrans.at(recv->idRecv()))
              for(GnssType typeObs : typesRecv)
                if(!typeObs.isInList(types.at(idAnt)))
                {
                  UInt idx;
                  if(typeObs.isInList(typesPattern, idx))
                    types.at(idAnt).push_back(typesPattern.at(idx));
                  else if(addNonMatchingTypes)
                    types.at(idAnt).push_back(typeObs & ~GnssType::PRN);
                }
          } // for(antenna)
      } // for(recv)

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
      logInfo<<countPara%"%9i receiver antenna parameters"s<<Log::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationReceiverAntennas::designMatrix(const GnssNormalEquationInfo &/*normalEquationInfo*/, const GnssObservationEquation &eqn, GnssDesignMatrix &A) const
{
  try
  {
    const UInt idAnt = receiver2antenna.at(eqn.receiver->idRecv());
    if((idAnt >= index.size()) || !index.at(idAnt).size())
      return;

    Matrix B(eqn.l.rows(), types.at(idAnt).size());
    std::vector<Bool> used(types.at(idAnt).size(), FALSE);
    UInt idPattern;
    for(UInt idType=0; idType<eqn.types.size(); idType++)
      if(eqn.types.at(idType).isInList(types.at(idAnt), idPattern))
      {
        axpy(1., eqn.A.column(GnssObservationEquation::idxUnit + idType), B.column(idPattern));
        used.at(idPattern) = TRUE;
      }
    if(!std::any_of(used.begin(), used.end(), [](Bool x){return x;}))
      return;

    const Matrix Acv = parametrization->designMatrix(eqn.azimutRecvAnt, eqn.elevationRecvAnt);
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
