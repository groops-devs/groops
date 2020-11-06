/***********************************************/
/**
* @file gnssParametrizationAntenna.cpp
*
* @brief GNSS antenna center variations.
*
* @author Torsten Mayer-Guerr
* @date 2019-05-31
*
*/
/***********************************************/

#define DOCSTRING_GnssParametrizationAntenna

#include "base/import.h"
#include "config/configRegister.h"
#include "files/fileGnssStationInfo.h"
#include "classes/parametrizationGnssAntenna/parametrizationGnssAntenna.h"
#include "gnss/gnss.h"
#include "gnss/gnssObservation.h"
#include "gnss/gnssDesignMatrix.h"
#include "gnss/gnssReceiver.h"
#include "gnss/gnssTransmitter.h"
#include "gnss/gnssParametrizationAntenna.h"

/***********************************************/

GROOPS_REGISTER_CLASS_WITHOUT_SUBS(GnssParametrizationAntenna, "gnssParametrizationAntennaType")
GROOPS_READCONFIG_CLASS(GnssParametrizationAntenna, "gnssParametrizationAntennaType")

/***********************************************/

GnssParametrizationAntenna::GnssParametrizationAntenna(Config &config, const std::string &name)
{
  try
  {
    readConfigSequence(config, name, Config::MUSTSET, "", "");
    readConfig(config, "antennaCenterVariations", parametrization,     Config::MUSTSET,  "",  "estimate antenna center variations");
    readConfig(config, "patternTypes",            typesPattern,        Config::OPTIONAL, "",  "gnssType for each pattern (first match is used)");
    readConfig(config, "addNonMatchingTypes",     addNonMatchingTypes, Config::DEFAULT,  "1", "add patterns for additional observed gnssTypes that don't match any of the above");
    readConfig(config, "groupAntennas",           ignoreSerial,        Config::DEFAULT,  "0", "common ACVs for same antenna build types (ignores antenna serial number)");
    endSequence(config);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationAntenna::initIntervalEarly(Gnss::AnalysisType /*analysisType*/, const std::vector<Time> &/*times*/, const Time &/*timeMargin*/, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    device2antenna.clear();
    antennaNames.clear();
    indexParameter.clear();
    types.clear();
    usedTypes.clear();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationAntenna::addAntennas(UInt idDevice, Bool deviceIsReceiver, const GnssStationInfo &stationInfo, const std::vector<GnssType> &typesObs)
{
  try
  {
    this->deviceIsReceiver = deviceIsReceiver;
    device2antenna.resize(deviceIsReceiver ? gnss().receiver.size() : gnss().transmitter.size(), NULLINDEX);

    Bool found = FALSE;
    for(const auto &antenna : stationInfo.antenna)
      if((antenna.timeEnd > gnss().times.at(0)) && (antenna.timeStart <= gnss().times.back()))
      {
        if(found)
        {
          logWarning<<gnss().times.at(0).dateTimeStr()<<" "<<stationInfo.markerName<<" has multiple antennas in interval -> only first one is used for parametrization/estimation"<<Log::endl;
          continue;
        }
        found = TRUE;

        // create name
        const std::string name = GnssAntennaDefinition::str(antenna.name, ignoreSerial ? ""s : antenna.serial, antenna.radome);

        // not already in list?
        const UInt idAnt = std::distance(antennaNames.begin(), std::find(antennaNames.begin(), antennaNames.end(), name));
        if(idAnt >= antennaNames.size())
        {
          antennaNames.push_back(name);
          types.resize(antennaNames.size(), typesPattern);
          usedTypes.resize(antennaNames.size(), std::vector<Bool>(typesPattern.size(), FALSE));
        }
        device2antenna.at(idDevice) = idAnt;

        // add observed types
        for(GnssType typeObs : typesObs)
        {
          typeObs &= ~GnssType::PRN;
          const UInt idx = GnssType::index(types.at(idAnt), typeObs);
          if(idx != NULLINDEX)
            usedTypes.at(idAnt).at(idx) = TRUE;
          else if(addNonMatchingTypes)
          {
            types.at(idAnt).push_back(typeObs);
            usedTypes.at(idAnt).push_back(TRUE);
          }
        } // for(typeObs)
      } // for(antenna)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationAntenna::initParameter(Gnss::NormalEquationInfo &normalEquationInfo)
{
  try
  {
    indexParameter.clear();

    if(!antennaNames.size())
      return;
    if(( deviceIsReceiver && (normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_RECEIVER_ANTENNACENTERVARIATIONS)) ||
       (!deviceIsReceiver && (normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_TRANSMITTER_ANTENNACENTERVARIATIONS)))
    {
      indexParameter.resize(antennaNames.size());
      for(UInt idAnt=0; idAnt<indexParameter.size(); idAnt++)
      {
        // remove unused types
        for(UInt i=types.at(idAnt).size(); i-->0;)
          if(!usedTypes.at(idAnt).at(i))
          {
            types.at(idAnt).erase(types.at(idAnt).begin()+i);
            usedTypes.at(idAnt).erase(usedTypes.at(idAnt).begin()+i);
          }

        // parameter names
        if(types.at(idAnt).size())
        {
          std::vector<ParameterName> baseNames;
          parametrization->parameterName(baseNames);
          for(GnssType type : types.at(idAnt))
          {
            std::string typeStr = "." + type.str();
            std::vector<ParameterName> parameterNames;
            for(const auto &base : baseNames)
              parameterNames.push_back( ParameterName(antennaNames.at(idAnt), base.type + typeStr, base.temporal, base.interval) );
            indexParameter.at(idAnt).push_back(normalEquationInfo.parameterNamesOther(parameterNames));
          }
        }
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool GnssParametrizationAntenna::isDesignMatrix(const Gnss::NormalEquationInfo &/*normalEquationInfo*/, UInt idRecv, UInt idTrans, UInt /*idEpoch*/) const
{
  try
  {
    const UInt idAnt = device2antenna.at(deviceIsReceiver ? idRecv : idTrans);
    return (idAnt < indexParameter.size()) && indexParameter.at(idAnt).size();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationAntenna::designMatrix(const Gnss::NormalEquationInfo &/*normalEquationInfo*/, const Gnss::ObservationEquation &eqn, Gnss::DesignMatrix &A) const
{
  try
  {
    if(!indexParameter.size())
      return;

    const UInt idAnt = device2antenna.at(deviceIsReceiver ? eqn.receiver->idRecv() : eqn.transmitter->idTrans());
    if((idAnt >= indexParameter.size()) || !indexParameter.at(idAnt).size())
      return;

    Matrix Acv;
    Matrix B(eqn.l.rows(), types.at(idAnt).size());
    std::vector<Bool> used(types.at(idAnt).size(), FALSE);
    if(deviceIsReceiver)
    {
      for(UInt idType=0; idType<eqn.types.size(); idType++)
        for(UInt idPattern=0; idPattern<types.at(idAnt).size(); idPattern++)
          if(types.at(idAnt).at(idPattern) == eqn.types.at(idType))
          {
            axpy(1., eqn.A.column(Gnss::ObservationEquation::idxUnit + idType), B.column(idPattern));
            used.at(idPattern) = TRUE;
            break;
          }
      if(!std::any_of(used.begin(), used.end(), [](Bool x){return x;}))
        return;
      Acv = parametrization->designMatrix(eqn.azimutRecvAnt, eqn.elevationRecvAnt);
    }
    else
    {
      for(UInt idType=0; idType<eqn.typesTransmitted.size(); idType++)
        for(UInt idPattern=0; idPattern<types.at(idAnt).size(); idPattern++)
          if(types.at(idAnt).at(idPattern) == eqn.typesTransmitted.at(idType))
          {
            axpy(1., eqn.A.column(Gnss::ObservationEquation::idxUnit + eqn.types.size() + idType), B.column(idPattern));
            used.at(idPattern) = TRUE;
            break;
          }
      if(!std::any_of(used.begin(), used.end(), [](Bool x){return x;}))
        return;
      Acv = parametrization->designMatrix(eqn.azimutTrans, eqn.elevationTrans);
    }

    for(UInt idPattern=0; idPattern<types.at(idAnt).size(); idPattern++)
      if(used.at(idPattern))
        matMult(1., B.column(idPattern), Acv, A.column(indexParameter.at(idAnt).at(idPattern)));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
