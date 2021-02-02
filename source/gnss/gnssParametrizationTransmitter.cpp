/***********************************************/
/**
* @file gnssParametrizationTransmitter.cpp
*
* @brief transmitter satellite systems.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2018-08-07
*
*/
/***********************************************/

#define DOCSTRING_GnssParametrizationTransmitter

#include "base/import.h"
#include "base/polynomial.h"
#include "base/kepler.h"
#include "base/planets.h"
#include "parser/dataVariables.h"
#include "config/configRegister.h"
#include "config/config.h"
#include "inputOutput/logging.h"
#include "files/fileMatrix.h"
#include "files/fileStringTable.h"
#include "files/fileGnssStationInfo.h"
#include "files/fileGnssSignalBias.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/parametrizationGravity/parametrizationGravity.h"
#include "classes/parametrizationAcceleration/parametrizationAcceleration.h"
#include "classes/timeSeries/timeSeries.h"
#include "gnss/gnss.h"
#include "gnss/gnssObservation.h"
#include "gnss/gnssDesignMatrix.h"
#include "gnss/gnssReceiver.h"
#include "gnss/gnssParametrizationEarthRotation.h"
#include "gnss/gnssParametrizationGravityField.h"
#include "gnss/gnssParametrizationTransmitter.h"
#include "gnss/gnssParametrizationTransmitterGps.h"
#include "gnss/gnssParametrizationTransmitterGlonass.h"
#include "gnss/gnssParametrizationTransmitterGalileo.h"
#include "gnss/gnssParametrizationTransmitterBeidou.h"
#include "gnss/gnssParametrizationTransmitter.h"

/***********************************************/

GROOPS_REGISTER_CLASS(GnssParametrizationTransmitter, "gnssParametrizationTransmitterType",
                      GnssParametrizationTransmitterGps,
                      GnssParametrizationTransmitterGlonass,
                      GnssParametrizationTransmitterGalileo,
                      GnssParametrizationTransmitterBeidou)

GROOPS_RENAMED_CLASS(gnssTransmitterConstellationType, gnssParametrizationTransmitterType, date2time(2020, 6, 12))

GROOPS_READCONFIG_CLASS(GnssParametrizationTransmitter, "gnssParametrizationTransmitterType")

/***********************************************/

GnssParametrizationTransmitterPtr GnssParametrizationTransmitter::create(Config &config, const std::string &name)
{
  try
  {
    GnssParametrizationTransmitterPtr system;
    std::string type;
    readConfigChoice(config, name, type, Config::MUSTSET, "", "");
    if(readConfigChoiceElement(config, "GPS", type, ""))
      system = GnssParametrizationTransmitterPtr(new GnssParametrizationTransmitterGps(config));
    if(readConfigChoiceElement(config, "GLONASS", type, ""))
      system = GnssParametrizationTransmitterPtr(new GnssParametrizationTransmitterGlonass(config));
    if(readConfigChoiceElement(config, "Galileo", type, ""))
      system = GnssParametrizationTransmitterPtr(new GnssParametrizationTransmitterGalileo(config));
    if(readConfigChoiceElement(config, "BeiDou", type, ""))
      system = GnssParametrizationTransmitterPtr(new GnssParametrizationTransmitterBeidou(config));
    endChoice(config);
    return system;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> Bool readConfig(Config &config, const std::string &name, GnssParametrizationTransmitter::BiasModel &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  try
  {
    if(!readConfigSequence(config, name, mustSet, defaultValue, annotation))
      return FALSE;

    renameDeprecatedConfig(config, "temporal", "parametrizationTemporal", date2time(2020, 8, 20));

    readConfig(config, "outputfileBias",          var.fileNameOut, Config::OPTIONAL, "",   "variable {prn} available");
    readConfig(config, "type",                    var.type,        Config::MUSTSET,  "",   "");
    readConfig(config, "parametrizationTemporal", var.temporal,    Config::MUSTSET,  "",   "");
    endSequence(config);
    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> Bool readConfig(Config &config, const std::string &name, GnssParametrizationTransmitter::TimeVariableBias &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  if(!readConfigSequence(config, name, mustSet, defaultValue, annotation))
    return FALSE;
  readConfig(config, "inputfileSignalBias",  var.inNameBias,  Config::MUSTSET, "", "variable {prn} available; instrument file with bias [m]");
  readConfig(config, "type",                 var.type,        Config::MUSTSET, "", "bias type");
  endSequence(config);
  return TRUE;
}

/***********************************************/
/***********************************************/

GnssParametrizationTransmitter::GnssParametrizationTransmitter()
{
  try
  {
    estimateCodeBias                = FALSE;
    estimatePhaseBias               = TRUE;
    estimateClockError              = EstimateClockError::NONE;
    estimateClockSigma              = FALSE;
    disableShadowEpochs             = FALSE;
    disablePostShadowRecoveryEpochs = FALSE;
    integrationDegree               = 7;

    estimateDynamicOrbit            = FALSE;
    minEstimableEpochsRatio         = 0.75;
    onlyEclipsingStochasticPulse    = FALSE;

    noPatternFoundAction = GnssAntennaDefinition::NoPatternFoundAction::USE_NEAREST_FREQUENCY;

    polynomial.init(7);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationTransmitter::init(const FileName &fileNameTransmitterList,
                                        const FileName &fileNameTransmitterInfo,
                                        const FileName &fileNameAntennaDef,
                                        const FileName &fileNameReceiverDef,
                                        UInt interpolationDegree)
{
  try
  {
    // read antenna defintion
    std::vector<GnssAntennaDefinitionPtr> antennaDefList;
    readFileGnssAntennaDefinition(fileNameAntennaDef, antennaDefList);

    // read receiver defintion
    std::vector<GnssReceiverDefinitionPtr> receiverDefList;
    if(!fileNameReceiverDef.empty())
      readFileGnssReceiverDefinition(fileNameReceiverDef, receiverDefList);

    // init transmitters
    // -----------------
    std::vector<std::string> transmitterList;
    readFileStringList(fileNameTransmitterList, transmitterList);

    VariableList fileNameVariableList;
    addVariable("prn", fileNameVariableList);

    for(UInt idTrans=0; idTrans<transmitterList.size(); idTrans++)
    {
      fileNameVariableList["prn"]->setValue(transmitterList.at(idTrans));

      GnssStationInfo transmitterInfo;
      readFileGnssStationInfo(fileNameTransmitterInfo(fileNameVariableList), transmitterInfo);
      transmitterInfo.fillAntennaPattern(antennaDefList);
      transmitterInfo.fillReceiverDefinition(receiverDefList);

      transmitter.push_back(std::make_shared<Transmitter>());
      transmitter.back()->base = this;
      transmitter.back()->transmitterInfo = transmitterInfo;
      transmitter.back()->type = GnssType("***"+transmitterInfo.markerNumber);
    }

    polynomial.init(interpolationDegree);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::vector<Gnss::TransmitterPtr> GnssParametrizationTransmitter::transmitters()
{
  std::vector<Gnss::TransmitterPtr> t;
  t.insert(t.begin(), transmitter.begin(), transmitter.end());
  return t;
}

/***********************************************/

std::vector<Gnss::ParametrizationPtr> GnssParametrizationTransmitter::parametrizations()
{
  std::vector<Gnss::ParametrizationPtr> p({shared_from_this()});
  p.insert(p.begin(), antennaCenterVariations.begin(), antennaCenterVariations.end());
  return p;
}

/***********************************************/

void GnssParametrizationTransmitter::initIntervalTransmitter(Gnss::AnalysisType /*analysisType*/, const std::vector<Time> &times, const Time &/*timeMargin*/, Parallel::CommunicatorPtr comm)
{
  try
  {
    this->times = times;

    // init transmitters
    // -----------------
    for(auto &trans : transmitter)
    {
      trans->useAtAll = TRUE;
      trans->use.resize(times.size(), TRUE);
      trans->signalBias = GnssSignalBias();
      trans->biasModel.clear();
    }

    VariableList fileNameVariableList;
    addVariable("prn", fileNameVariableList);

    // ===========================================================

    // init variational equations
    // --------------------------
    UInt countTrans = 0;
    for(auto &trans : transmitter)
    {
      try
      {
        fileNameVariableList["prn"]->setValue(trans->name());
        VariationalEquationFromFile file;

        // gravity field parametrization
        ParametrizationGravityPtr parametrizationGravity;
        if(gnss().gravityField)
          parametrizationGravity = gnss().gravityField->parametrization();

        // stochastic pulses in interval
        std::vector<Time> pulse = stochasticPulse;
        pulse.erase(std::remove_if(pulse.begin(), pulse.end(), [&](auto &p) {return (p <= times.at(0)) || (p >= times.back());}), pulse.end());

        // setup stochastic pulses only if satellite is eclipsing?
        if(pulse.size() && eclipse && onlyEclipsingStochasticPulse)
        {
          file.open(fileNameTemplateInVariational(fileNameVariableList), parametrizationGravity, parametrizationAcceleration, pulse, ephemerides, integrationDegree);
          auto variationalEquation = file.integrateArc(times.at(0), times.back(), TRUE/*computePosition*/, FALSE/*computeVelocity*/);
          Bool isEclipsing = FALSE;
          for(UInt idEpoch=0; idEpoch<variationalEquation.times.size(); idEpoch++)
            if(eclipse->factor(variationalEquation.times.at(idEpoch), Vector3d(variationalEquation.pos0.row(3*idEpoch, 3)), ephemerides) < 1)
            {
              isEclipsing = TRUE;
              break;
            }
          if(!isEclipsing)
            pulse.clear();
        }

        // read data from variational equations
        file.open(fileNameTemplateInVariational(fileNameVariableList), parametrizationGravity, parametrizationAcceleration, pulse, ephemerides, integrationDegree);
        trans->variationalEquation = file.integrateArc(times.at(0), times.back(), TRUE/*computePosition*/, TRUE/*computeVelocity*/);
        trans->satelliteModel      = file.satellite();

        // disable epochs that are outside variational orbit time period
        for(UInt idEpoch = 0; idEpoch < times.size(); idEpoch++)
          if(times.at(idEpoch) < trans->variationalEquation.times.at(0) ||
             times.at(idEpoch) > trans->variationalEquation.times.back())
            trans->disable(idEpoch);

        // parameter count
        trans->countParameterSatellite    = file.parameterCountSatellite();
        trans->countParameterSatelliteArc = file.parameterCountSatelliteArc();
        trans->xOrbit = Vector(trans->countParameterSatellite + trans->countParameterSatelliteArc);

        // parameter names
        trans->parameterNameSatellite.clear();
        file.parameterNameSatellite(trans->parameterNameSatellite);
        for(auto &name : trans->parameterNameSatellite)
          name.object = trans->name();

        trans->parameterNameSatelliteArc.clear();
        file.parameterNameSatelliteArc(trans->parameterNameSatelliteArc);
        for(auto &name : trans->parameterNameSatelliteArc)
          name.object = trans->name();

        countTrans++;
      }
      catch(std::exception &/*e*/)
      {
        trans->disable();
        continue;
      }
    }

    if(!countTrans && Parallel::isMaster(comm))
    {
      fileNameVariableList["prn"]->setValue("***");
      logWarning<<"Initialization of all variational equations failed. Wrong file name <"<<fileNameTemplateInVariational(fileNameVariableList)<<">?"<<Log::endl;
    }


    // interpolate satellite orientation from variationalEquation.times to times
    // -------------------------------------------------------------------------
    for(auto &trans : transmitter)
      if(trans->useable())
      {
        Matrix quaternion(trans->variationalEquation.times.size(), 4);
        for(UInt i=0; i<quaternion.rows(); i++)
          copy(trans->variationalEquation.rotSat.at(i).quaternion().trans(), quaternion.row(i));
        for(UInt i=1; i<quaternion.rows(); i++)
          if(inner(quaternion.row(i), quaternion.row(i-1))<0)
            quaternion.row(i) *= -1;
        quaternion = polynomial.interpolate(times, trans->variationalEquation.times, quaternion);
        for(UInt idEpoch=0; idEpoch<quaternion.rows(); idEpoch++)
          quaternion.row(idEpoch) *= 1./norm(quaternion.row(idEpoch));

        trans->crf2arf.resize(times.size());
        for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
        {
          const UInt idAntenna = trans->transmitterInfo.findAntenna(times.at(idEpoch));
          if(idAntenna >= trans->transmitterInfo.antenna.size())
            throw(Exception(trans->name() + ": no antenna found at " + times.at(idEpoch).dateTimeStr()));
          trans->crf2arf.at(idEpoch) = trans->transmitterInfo.antenna.at(idAntenna).local2antennaFrame * inverse(Rotary3d(quaternion.row(idEpoch).trans()));
        }
      }

    // exclusion of shadow crossing and post-shadow recovery epochs
    // ------------------------------------------------------------
    if(eclipse && (disableShadowEpochs || disablePostShadowRecoveryEpochs))
    {
      for(auto &trans : transmitter)
        if(trans->useable())
        {
          // 30 minute post-shadow recovery time for block IIA satellites
          Time recoveryTime;
          if(trans->transmitterInfo.antenna.at(trans->transmitterInfo.findAntenna(times.at(0))).name.find(std::string("IIA"))!=std::string::npos)
            recoveryTime = seconds2time(30*60);

          Time timeShadowExit;

          // Look for a potential shadow exit before the start of the interval that would result in the
          // satellite already performing a post-shadow recovery maneuver at the start of the interval
          if(disablePostShadowRecoveryEpochs && recoveryTime.seconds() > 0)
          {
            Kepler kepler(times.at(0), trans->positionCoM(0, times.at(0)), trans->velocity(0, times.at(0)));
            for(UInt i=0; i<recoveryTime.seconds()/60; i++)
            {
              const Time time = times.at(0) - seconds2time(i*60.);
              if(eclipse->factor(time, kepler.position(time), ephemerides) < 0.5)
              {
                timeShadowExit = time;
                break;
              }
            }
          }

          Double factorPreviousEpoch = 1.0;
          for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
          {
            const Double factor = eclipse->factor(times.at(idEpoch), trans->positionCoM(idEpoch, times.at(idEpoch)), ephemerides);
            if(factorPreviousEpoch < 0.5 && factor >= 0.5)
              timeShadowExit = times.at(idEpoch);

            // set satellite unuseable during shadow crossing and post-shadow recovery maneuver
            if((disableShadowEpochs && factor < 0.5) || (disablePostShadowRecoveryEpochs && times.at(idEpoch) < timeShadowExit+recoveryTime))
              trans->disable(idEpoch);

            factorPreviousEpoch = factor;
          } // for(idEpoch)
        }
    }

    // ===========================================================

    // satellite clocks
    // ----------------
    for(auto &trans : transmitter)
      if(trans->useable())
      {
        MiscValueArc arc;
        fileNameVariableList["prn"]->setValue(trans->name());
        try
        {
          arc = InstrumentFile::read(fileNameTemplateInClock(fileNameVariableList));
          if(!arc.size())
            throw(Exception("empty file"));
        }
        catch(std::exception &/*e*/)
        {
          if(Parallel::isMaster(comm))
            logWarning<<"Unable to read clock file <"<<fileNameTemplateInClock(fileNameVariableList)<<">, disabling transmitter."<<Log::endl;
          trans->disable();
          continue;
        }

        trans->clock0 = Vector(times.size());
        trans->dclock = Vector(times.size());

        std::vector<Time> timesArc = arc.times();
        const Double clockSampling = medianSampling(timesArc).seconds();
        UInt idx = 0;
        for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
        {
          // find interval
          while((idx+1<timesArc.size()) && (timesArc.at(idx+1)<=times.at(idEpoch)))
            idx++;

          if(timesArc.at(idx) == times.at(idEpoch))
            trans->dclock(idEpoch) = arc.at(idx).value * LIGHT_VELOCITY;
          else if(idx+1<timesArc.size() && ((estimateClockError != EstimateClockError::NONE) || (timesArc.at(idx+1)-timesArc.at(idx)).seconds()<1.5*clockSampling))
          {
            // linear interpolation
            const Double t = (times.at(idEpoch)-timesArc.at(idx)).seconds()/(timesArc.at(idx+1)-timesArc.at(idx)).seconds();
            trans->dclock(idEpoch) = ((1-t)*arc.at(idx).value + t*arc.at(idx+1).value)  * LIGHT_VELOCITY;
          }
          else
            trans->disable(idEpoch);
        }
      }

    // ===========================================================

    // code & phase biases
    // -------------------
    if(!fileNameTemplateInSignal.empty())
      for(auto &trans : transmitter)
        if(trans->useable())
        {
          fileNameVariableList["prn"]->setValue(trans->name());
          try
          {
            readFileGnssSignalBias(fileNameTemplateInSignal(fileNameVariableList), trans->signalBias);
          }
          catch(std::exception &/*e*/)
          {
            if(Parallel::isMaster(comm))
              logWarning<<"Unable to read signal bias file <"<<fileNameTemplateInSignal(fileNameVariableList)<<">, disabling transmitter."<<Log::endl;
            trans->disable();
          }
        }

    // ===========================================================

    // time-variable signal biases
    // ---------------------------
    for(const auto &bias : timeVariableBiases)
    {
      Bool foundAnyFile = FALSE;
      for(auto &trans : transmitter)
        if(trans->useable())
        {
          try
          {
            fileNameVariableList["prn"]->setValue(trans->name());
            MiscValueArc arc = InstrumentFile::read(bias.inNameBias(fileNameVariableList));
            TimeVariableBias tvBias;
            tvBias.type = bias.type;
            tvBias.biases = polynomial.interpolate(times, arc.times(), arc.matrix().column(1));
            trans->timeVariableBiases.push_back(tvBias);
            foundAnyFile = TRUE;
          }
          catch(std::exception &/*e*/)
          {
            //if(Parallel::isMaster(comm))
            //  logWarning<<trans->name()<<": Unable to read time-variable signal bias file for "<<bias.type.str()<<"."<<Log::endl;
          }
        }
      if(Parallel::isMaster(comm) && !foundAnyFile)
        logWarning<<name()<<": "<<bias.type.str()<<" time-variable signal bias file not found for any transmitter."<<Log::endl;
    }

    // ===========================================================

    if(!fileNameClockNoise.empty())
      for(auto &trans : transmitter)
        if(trans->useable())
        {
          trans->clockSigma = 1.;
          try
          {
            fileNameVariableList["prn"]->setValue(trans->name());
            readFileMatrix(fileNameClockNoise(fileNameVariableList), trans->clockNoise);
          }
          catch(std::exception &/*e*/)
          {
            if(Parallel::isMaster(comm))
              logWarning<<trans->name()<<": Unable to read clock noise file, disabling transmitter."<<Log::endl;
            trans->disable();
          }
        }

    // ===========================================================

    // test antennas
    // -------------
    for(const auto &trans : transmitter)
      if(trans->useable())
      {
        for(const auto &antenna : trans->transmitterInfo.antenna)
          if(Parallel::isMaster(comm) && antenna.timeEnd >= times.front() && antenna.timeStart < times.back() && !antenna.antennaDef)
            logWarning<<trans->name()<<": No antenna definition found for "<<antenna.name<<"."<<antenna.serial<<Log::endl;

        for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
        {
          UInt idAnt = trans->transmitterInfo.findAntenna(times.at(idEpoch));
          if(idAnt == NULLINDEX || !trans->transmitterInfo.antenna.at(idAnt).antennaDef)
            trans->disable(idEpoch);
        }
      }

    // test useable
    // ------------
    for(auto &trans : transmitter)
      if(trans->useable())
      {
        Bool useAtAll = FALSE;
        for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
          if(trans->useable(idEpoch))
            useAtAll = TRUE;
        if(!useAtAll)
        {
          if(Parallel::isMaster(comm))
            logWarning<<trans->name()<<": No usable epochs, disabling transmitter."<<Log::endl;
          trans->disable();
        }
      }

    // ===========================================================

    // init temporal parametrization
    // -----------------------------
    if(clockModel)
    {
      clockModel->setInterval(times.front(), times.back(), TRUE);
      // move reference clocks to clock0
      for(auto &trans : transmitter)
        if(trans->useable())
        {
          trans->clock0 = trans->dclock;
          Matrix A(times.size(), clockModel->parameterCount());
          for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
            if(trans->useable(idEpoch))
              copy(clockModel->factors(times.at(idEpoch)).trans(), A.row(idEpoch));
          reduceLeastSquaresFit(A, trans->dclock);
          trans->clock0 -= trans->dclock;
        }
    }

    for(auto &model : biasModel)
    {
      model.x    = Vector(model.temporal->parameterCount());
      model.bias = Vector(times.size());
      model.indexParameter = Gnss::ParameterIndex();
      model.temporal->setInterval(times.front(), times.back(), TRUE);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool GnssParametrizationTransmitter::initIntervalDisabling(Gnss::AnalysisType /*analysisType*/, const std::vector<Time> &/*times*/, const Time &/*timeMargin*/, Parallel::CommunicatorPtr comm)
{
  try
  {
    // check available observations
    // ----------------------------
    Bool disabled = FALSE;
    for(auto &trans : transmitter)
      if(trans->useable())
      {
        UInt epochCount = 0;
        for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
          if(trans->useable(idEpoch))
          {
            UInt count = 0;
            for(const auto &recv : gnss().receiver)
              if(recv->useable())
                count += recv->countObservations(trans->idTrans(), idEpoch, idEpoch);
            Parallel::reduceSum(count, 0, comm);
            Parallel::broadCast(count, 0, comm);
            if(count == 0)
            {
              trans->disable(idEpoch);
              disabled = TRUE;
            }
            else
              epochCount++;
          }

        if(estimateDynamicOrbit && (epochCount < minEstimableEpochsRatio*times.size()))
        {
          if(Parallel::isMaster(comm))
            logWarning<<trans->name()<<" disabled: not enough estimable epochs ("<<epochCount<<") to estimate dynamic orbit"<<Log::endl;
          trans->disable();
          disabled = TRUE;
        }
      }

    return disabled;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationTransmitter::initIntervalLate(Gnss::AnalysisType /*analysisType*/, const std::vector<Time> &/*times*/, const Time &/*timeMargin*/, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    VariableList fileNameVariableList;
    addVariable("prn", fileNameVariableList);

    for(auto &trans : transmitter)
      if(trans->useable())
      {
        // check available observations types
        // ----------------------------------
        std::vector<GnssType> types;
        for(auto &typesTrans : gnss().typesRecvTrans)
          for(GnssType type : typesTrans.at(trans->idTrans()))
            if((type == GnssType::PHASE) || (type == GnssType::RANGE))
              if(GnssType::index(types, type) == NULLINDEX)
                types.push_back(type + trans->PRN());
        types = gnss().replaceCompositeSignals(types);

        // code & phase biases
        // -------------------
        if(types.size()) // consider simulation case
        {
          trans->signalBias.bias = trans->signalBias.compute(types); // apriori signal bias
          trans->signalBias.type = types;
        }

        // signal bias model
        // -----------------
        fileNameVariableList["prn"]->setValue(trans->name());
        for(auto &model : biasModel)
          if(GnssType::index(types, model.type) != NULLINDEX)
            trans->biasModel.push_back(model);

        // antenna center variations
        // -------------------------
        for(auto &acv : antennaCenterVariations)
          acv->addAntennas(trans->idTrans(), FALSE/*deviceIsReceiver*/, trans->transmitterInfo, types);
      } // for(idTrans)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationTransmitter::initParameter(Gnss::NormalEquationInfo &normalEquationInfo)
{
  try
  {
    for(auto &trans : transmitter)
    {
      trans->indexParameterClock.clear();
      trans->indexParameterClock.resize(times.size(), Gnss::ParameterIndex());
      trans->indexParameterClockModel   = Gnss::ParameterIndex();
      trans->indexParameterSatellite    = Gnss::ParameterIndex();
      trans->indexParameterSatelliteArc = Gnss::ParameterIndex();
      trans->indexParameterBias         = Gnss::ParameterIndex();
      for(auto &model : trans->biasModel)
        model.indexParameter = Gnss::ParameterIndex();
    }

    // transmitter
    // -----------
    for(auto trans : transmitter)
      if(trans->useable())
      {
        // epoch parameter
        // ---------------
        for(UInt idEpoch : normalEquationInfo.idEpochs)
          if(trans->useable(idEpoch))
          {
            // clock
            if((estimateClockError != EstimateClockError::NONE) && (normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_TRANSMITTER_CLOCK))
            {
              // check if clock parameter is estimable with currently selected receivers
              UInt count = 0;
              for(UInt idRecv = 0; idRecv < gnss().receiver.size(); idRecv++)
                if(gnss().receiver.at(idRecv)->useable() && normalEquationInfo.estimateReceiver.at(idRecv))
                  count += gnss().receiver.at(idRecv)->countObservations(trans->idTrans(), idEpoch, idEpoch);
              Parallel::reduceSum(count, 0, normalEquationInfo.comm);
              Parallel::broadCast(count, 0, normalEquationInfo.comm);
              if(!count)
                continue;

              trans->indexParameterClock.at(idEpoch) = normalEquationInfo.parameterNamesEpochTransmitter(idEpoch, trans->idTrans(), {ParameterName(trans->name(), "clock", "", times.at(idEpoch))});
            }
          }

        // variational equations
        // ---------------------
        if(estimateDynamicOrbit && (normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_TRANSMITTER_ORBIT))
        {
          trans->indexParameterSatellite    = normalEquationInfo.parameterNamesTransmitter(trans->idTrans(), trans->parameterNameSatellite);
          trans->indexParameterSatelliteArc = normalEquationInfo.parameterNamesTransmitter(trans->idTrans(), trans->parameterNameSatelliteArc);
        } // if(estimateDynamicOrbit)

        // clock model
        // -----------
        if(clockModel && clockModel->parameterCount() && (normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_TRANSMITTER_CLOCKMODEL))
        {
          std::vector<ParameterName> parameterNames;
          clockModel->parameterName({ParameterName(trans->name(), "clock")}, parameterNames);
          trans->indexParameterClockModel = normalEquationInfo.parameterNamesTransmitter(trans->idTrans(), parameterNames);
        }

        // signal bias
        // -----------
        if((normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_TRANSMITTER_SIGNALBIAS) ||
           (normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_TRANSMITTER_TECBIAS))
        {
          // determine groups: for each receiver: used types (receiver types)
          std::vector<std::vector<std::vector<GnssType>>> typesRecvTrans(gnss().receiver.size(), std::vector<std::vector<GnssType>>(1));
          for(UInt idRecv=0; idRecv<gnss().receiver.size(); idRecv++)
            if(normalEquationInfo.estimateReceiver.at(idRecv))
              typesRecvTrans.at(idRecv).at(0) = gnss().typesRecvTrans.at(idRecv).at(trans->idTrans());
          auto groupTypes = gnss().groupTypes(typesRecvTrans, 0, FALSE/*isReceiver*/);

          if(groupTypes.size())
          {
            std::vector<ParameterName> parameterNames;
            gnss().signalBiasParameter(trans->signalBias.type, gnss().groupTypes(typesRecvTrans, 0, FALSE/*isReceiver*/),
                                      /*eliminateClock*/(normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_TRANSMITTER_CLOCK),
                                      trans->isCodeBiasEstimated(normalEquationInfo),
                                      trans->isPhaseBiasEstimated(normalEquationInfo),
                                      /*estimateTecBias*/(normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_TRANSMITTER_TECBIAS),
                                      trans->Bias, parameterNames);
            if(parameterNames.size())
            {
              for(auto &parameterName : parameterNames)
                parameterName.object = trans->name();
              trans->indexParameterBias = normalEquationInfo.parameterNamesTransmitter(trans->idTrans(), parameterNames);
            }
          }
        }

        // bias model
        // -----------
        if(normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_TRANSMITTER_SIGNALBIASMODEL)
        {
          for(auto &model : trans->biasModel)
            if(model.temporal && model.temporal->parameterCount())
            {
              std::vector<ParameterName> parameterNames;
              model.temporal->parameterName({ParameterName(trans->name(), "signalBias."+model.type.str())}, parameterNames);
              model.indexParameter = normalEquationInfo.parameterNamesTransmitter(trans->idTrans(), parameterNames);
            }
        }
      } // for(trans)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationTransmitter::aprioriParameter(const Gnss::NormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const
{
  try
  {
    if(!Parallel::isMaster(normalEquationInfo.comm))
      return;

    for(auto &trans : transmitter)
      if(trans->useable())
      {
        // clock error
        // -----------
        for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
          if(trans->indexParameterClock.at(idEpoch))
            x0(normalEquationInfo.index(trans->indexParameterClock.at(idEpoch)), 0) = LIGHT_VELOCITY * trans->clockError(idEpoch, times.at(idEpoch));

        // update positions & velocities
        // -----------------------------
        if(trans->indexParameterSatellite)
          copy(trans->xOrbit.row(0, trans->countParameterSatellite),
               x0.row(normalEquationInfo.index(trans->indexParameterSatellite), trans->countParameterSatellite));
        if(trans->indexParameterSatelliteArc)
          copy(trans->xOrbit.row(trans->countParameterSatellite, trans->countParameterSatelliteArc),
               x0.row(normalEquationInfo.index(trans->indexParameterSatelliteArc), trans->countParameterSatelliteArc));

        // signal biases
        // -------------
        if(trans->indexParameterBias)
        {
          Vector l = trans->signalBias.bias;
          Matrix A = trans->Bias;
          copy(leastSquares(A, l), x0.row(normalEquationInfo.index(trans->indexParameterBias), trans->Bias.columns()));
        }

        // signal bias model
        // -----------------
        for(auto &model : trans->biasModel)
          if(model.indexParameter)
          {
            Matrix A(times.size(), model.temporal->parameterCount());
            for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
              copy(model.temporal->factors(times.at(idEpoch)).trans(), A.row(idEpoch));
            copy(leastSquares(A, Vector(model.bias)), x0.row(normalEquationInfo.index(model.indexParameter), model.temporal->parameterCount()));
          } // for(model)
      } // for(trans)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationTransmitter::observationEquationEpoch(const Gnss::NormalEquationInfo &normalEquationInfo, UInt idEpoch, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const
{
  try
  {
    // zero-mean constraint of satellite clocks
    // ----------------------------------------
    if(Parallel::isMaster(normalEquationInfo.comm) &&
       sigmaClockZeroMean && (estimateClockError == EstimateClockError::EPOCH) &&
       (normalEquationInfo.estimationType & Gnss::NormalEquationInfo::CONSTRAINT_TRANSMITTER_CLOCK_ZEROMEAN))
    {
      const UInt count = std::count_if(transmitter.begin(), transmitter.end(), [idEpoch](auto &trans) {return trans->indexParameterClock.at(idEpoch);});
      Gnss::DesignMatrix A(normalEquationInfo, Vector(1));
      for(auto &trans : transmitter)
        if(trans->indexParameterClock.at(idEpoch))
        {
//           A.l(0) -= trans->dclock(idEpoch)/count/sigma; // remove apriori value -> regularization towards 0
          A.column(trans->indexParameterClock.at(idEpoch))(0,0) = 1.0/count/sigmaClockZeroMean;
        }
      A.accumulateNormals(normals, n, lPl, obsCount);
    } // if(estimateClockError)

    // clock noise model
    // -----------------
    if(Parallel::isMaster(normalEquationInfo.comm) && (estimateClockError == EstimateClockError::NOISEMODEL) &&
       (normalEquationInfo.estimationType & Gnss::NormalEquationInfo::CONSTRAINT_TRANSMITTER_CLOCKMODEL))
    {
      const UInt countTrans = std::count_if(transmitter.begin(), transmitter.end(), [idEpoch](auto &trans) {return trans->indexParameterClock.at(idEpoch);});
      if(countTrans)
      {
        Gnss::DesignMatrix A(normalEquationInfo, Vector(countTrans));
        UInt idTrans = 0;
        for(auto &trans : transmitter)
          if(trans->indexParameterClock.at(idEpoch))
          {
            // how many epochs back in time?
            UInt count = 0;
            while((count < trans->clockNoise.rows()) && (count <= idEpoch) && (trans->indexParameterClock.at(idEpoch-count)))
              count++;

            for(UInt i=0; i<count; i++)
            {
              A.l(idTrans) -= trans->clockNoise(count-1, count-1-i) * trans->dclock(idEpoch-i); // remove apriori value -> regularization towards 0
              A.column(trans->indexParameterClock.at(idEpoch-i))(idTrans,0) = trans->clockNoise(count-1, count-1-i);
            }
            idTrans++;
          }
        A.accumulateNormals(normals, n, lPl, obsCount);
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationTransmitter::observationEquation(const Gnss::NormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const
{
  try
  {
    if(!Parallel::isMaster(normalEquationInfo.comm))
      return;

    // ============================================================

    // determine null space of signal biases
    // -------------------------------------
    if((normalEquationInfo.estimationType & Gnss::NormalEquationInfo::CONSTRAINT_TRANSMITTER_SIGNALBIAS_ZEROMEAN) &&
       std::any_of(transmitter.begin(), transmitter.end(), [](const auto &trans){return trans->indexParameterBias;}))
    {
      // simulate simplified normal equations (each receiver group observes all satellites at onces)
      // parameters: all transmitter biases
      Matrix N;
      std::vector<GnssType> types;
      {
        // for each receiver and transmitter: used types (receiver types)
        std::vector<std::vector<std::vector<GnssType>>> typesRecvTransAll(gnss().receiver.size(), std::vector<std::vector<GnssType>>(transmitter.size()));
        for(UInt idRecv=0; idRecv<gnss().receiver.size(); idRecv++)
          if(normalEquationInfo.estimateReceiver.at(idRecv))
            for(UInt idTrans=0; idTrans<transmitter.size(); idTrans++)
              typesRecvTransAll.at(idRecv).at(idTrans) = gnss().typesRecvTrans.at(idRecv).at(transmitter.at(idTrans)->idTrans());
        // keep only phase and range
        for(auto &typesRecv : typesRecvTransAll)
          for(auto &types : typesRecv)
            types.erase(std::remove_if(types.begin(), types.end(), [](GnssType &t) {return (t != GnssType::PHASE) && (t != GnssType::RANGE);}), types.end());
        // summarize similar receivers
        std::vector<std::vector<std::vector<GnssType>>> typesRecvTrans;
        for(UInt idRecv=0; idRecv<gnss().receiver.size(); idRecv++)
          if(std::any_of(typesRecvTransAll.at(idRecv).begin(), typesRecvTransAll.at(idRecv).end(), [](const auto &group){return group.size();}))
            if(std::find(typesRecvTrans.begin(), typesRecvTrans.end(), typesRecvTransAll.at(idRecv)) == typesRecvTrans.end())
              typesRecvTrans.push_back(typesRecvTransAll.at(idRecv));
        // count available observations types
        for(auto &typesRecv : typesRecvTrans)
          for(auto &typesTrans : typesRecv)
            for(GnssType type : typesTrans)
              if(GnssType::index(types, type) == NULLINDEX)
                types.push_back(type);
        types = gnss().replaceCompositeSignals(types);

        // transmitter biases
        UInt countTrans = 0;
        std::vector<UInt> idxBiasTrans(1, 0);
        for(auto &trans : transmitter)
          if(trans->indexParameterBias)
          {
            idxBiasTrans.push_back(idxBiasTrans.back() + trans->Bias.columns());
            countTrans++;
          }

        // receiver biases
        std::vector<std::vector<GnssType>> typesRecv(typesRecvTrans.size());
        std::vector<Matrix>                BiasRecvList(typesRecvTrans.size());
        std::vector<UInt>                  idxBiasRecv(1, 0);
        for(UInt idRecv=0; idRecv<typesRecvTrans.size(); idRecv++)
        {
          for(auto &typesTrans : typesRecvTrans.at(idRecv))
            for(GnssType type : typesTrans)
              if(GnssType::index(typesRecv.at(idRecv), type) == NULLINDEX)
                typesRecv.at(idRecv).push_back(type);
          std::vector<ParameterName> parameterNamesRecv;
          gnss().signalBiasParameter(typesRecv.at(idRecv), gnss().groupTypes(typesRecvTrans, idRecv, TRUE/*isReceiver*/),
                                     /*eliminateClock*/   (normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_RECEIVER_CLOCK),
                                     /*estimateCodeBias*/ (normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_RECEIVER_SIGNALBIAS),
                                     /*estimatePhaseBias*/(normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_RECEIVER_SIGNALBIAS),
                                     /*estimateTecBias*/  (normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_RECEIVER_TECBIAS),
                                     BiasRecvList.at(idRecv), parameterNamesRecv);
          idxBiasRecv.push_back(idxBiasRecv.back() + BiasRecvList.at(idRecv).columns());
        }

        // receiver clocks
        const UInt countClockRecv = (normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_RECEIVER_CLOCK) ? typesRecvTrans.size() : 0;
        const UInt idxClockRecv   = idxBiasRecv.back();

        // transmitter clocks
        const UInt countClockTrans = (normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_TRANSMITTER_CLOCK) ? countTrans : 0;
        const UInt idxClockTrans   = idxBiasRecv.back() + countClockRecv;

        // accumulate normal equations
        Matrix N11(idxBiasRecv.back() + countClockRecv + countClockTrans, Matrix::SYMMETRIC); // receiver biases and receiver/transmitter clocks
        Matrix N12(idxBiasRecv.back() + countClockRecv + countClockTrans, idxBiasTrans.back());
        Matrix N22(idxBiasTrans.back(), Matrix::SYMMETRIC); // transmitter biases
        UInt idxTrans = 0;
        for(UInt idTrans=0; idTrans<transmitter.size(); idTrans++)
          if(transmitter.at(idTrans)->indexParameterBias)
          {
            for(UInt idRecv=0; idRecv<typesRecvTrans.size(); idRecv++)
              if(typesRecvTrans.at(idRecv).at(idTrans).size())
              {
                const std::vector<GnssType> &types = typesRecvTrans.at(idRecv).at(idTrans);

                // receiver biases
                Matrix BiasRecv(types.size(), BiasRecvList.at(idRecv).columns());
                for(UInt idType=0; idType<types.size(); idType++)
                  copy(BiasRecvList.at(idRecv).row(GnssType::index(typesRecv.at(idRecv), types.at(idType))), BiasRecv.row(idType));

                // transmitter biases
                Matrix BiasTrans(types.size(), transmitter.at(idTrans)->Bias.columns());
                Matrix T;
                std::vector<GnssType> typesTrans;
                gnss().defaultSignalComposition(types, typesTrans, T);
                for(UInt idType=0; idType<typesTrans.size(); idType++)
                {
                  const UInt idx = GnssType::index(transmitter.at(idTrans)->signalBias.type, typesTrans.at(idType));
                  matMult(1., T.column(idType), transmitter.at(idTrans)->Bias.row(idx), BiasTrans);
                }

                // clocks
                Vector Clock(types.size(), 1.);

                // eliminate STEC
                if((normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_IONOSPHERE_STEC) &&
                  !(normalEquationInfo.estimationType & Gnss::NormalEquationInfo::CONSTRAINT_IONOSPHERE_STEC))
                {
                  Vector STEC(types.size());
                  for(UInt idType=0; idType<types.size(); idType++)
                    STEC(idType) = types.at(idType).ionosphericFactor();
                  eliminationParameter(STEC, {BiasTrans, BiasRecv, Clock});
                }

                // accumulate normals
                rankKUpdate(1, BiasRecv,                 N11.slice(idxBiasRecv.at(idRecv),    idxBiasRecv.at(idRecv),    BiasRecv.columns(),  BiasRecv.columns()));
                rankKUpdate(1, BiasTrans,                N22.slice(idxBiasTrans.at(idxTrans), idxBiasTrans.at(idxTrans), BiasTrans.columns(), BiasTrans.columns()));
                matMult(1., BiasRecv.trans(), BiasTrans, N12.slice(idxBiasRecv.at(idRecv),    idxBiasTrans.at(idxTrans), BiasRecv.columns(),  BiasTrans.columns()));
                if(normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_RECEIVER_CLOCK)
                {
                  rankKUpdate(1, Clock,                 N11.slice(idxClockRecv+idRecv,    idxClockRecv+idRecv,       1, 1));
                  matMult(1., BiasRecv.trans(), Clock,  N11.slice(idxBiasRecv.at(idRecv), idxClockRecv+idRecv,       BiasRecv.columns(), 1));
                  matMult(1., Clock.trans(), BiasTrans, N12.slice(idxClockRecv+idRecv,    idxBiasTrans.at(idxTrans), 1, BiasTrans.columns()));
                }
                if(normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_TRANSMITTER_CLOCK)
                {
                  rankKUpdate(1, Clock,                 N11.slice(idxClockTrans+idxTrans,  idxClockTrans+idxTrans,    1, 1));
                  matMult(1., BiasRecv.trans(), Clock,  N11.slice(idxBiasRecv.at(idRecv),  idxClockTrans+idxTrans,    BiasRecv.columns(), 1));
                  matMult(1., Clock.trans(), BiasTrans, N12.slice(idxClockTrans+idxTrans,  idxBiasTrans.at(idxTrans), 1, BiasTrans.columns()));
                }
                if((normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_RECEIVER_CLOCK) && (normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_TRANSMITTER_CLOCK))
                  matMult(1, Clock.trans(), Clock, N11.slice(idxClockRecv+idRecv, idxClockTrans+idxTrans, 1, 1));
              }
            idxTrans++;
          }
        // zero mean of transmitter clocks
        rankKUpdate(1, Vector(countClockTrans, 1.).trans(), N11.slice(idxBiasRecv.back(), idxBiasRecv.back(), countClockTrans, countClockTrans));

        // eliminate bias and clock parameters (with pseudo inverse)
        {
          Vector eigen = eigenValueDecomposition(N11, TRUE/*computeEigenVectors*/);
          for(UInt i=0; i<eigen.rows(); i++)
            N11.column(i) *= (eigen(i) < 1e-8) ? 0. : (1./std::sqrt(eigen(i)));
          N12 = N11.trans() * N12;
          rankKUpdate(-1, N12, N22);
        }
        // count zero eigen values => null space
        Vector eigen = eigenValueDecomposition(N22, TRUE/*computeEigenVectors*/);
        UInt count = 0;
        while((count<eigen.rows()) && (eigen(count) < 1e-8))
          count++;
        // null space defines the constraint equations
        N = N22.column(0, count).trans();
      }

      logStatus<<"apply "<<N.rows()<<" zero mean equations for "<<types.size()<<" signal bias types of "<<name()<<" transmitters"<<Log::endl;
      if(N.rows())
      {
        const Double sigma = 0.000001; // 0.1 mm
        Gnss::DesignMatrix A(normalEquationInfo, Vector(N.rows()));
        UInt idx = 0;
        for(auto &trans : transmitter)
          if(trans->indexParameterBias)
          {
            axpy(1./sigma, N.column(idx, trans->Bias.columns()), A.column(trans->indexParameterBias));
            idx += trans->Bias.columns();
          }
        A.accumulateNormals(normals, n, lPl, obsCount);
      }
    }

    // ============================================================

    // temporal zero-mean constraint of signal bias models
    // ---------------------------------------------------
    if(normalEquationInfo.estimationType & Gnss::NormalEquationInfo::CONSTRAINT_TRANSMITTER_SIGNALBIASMODEL_ZEROMEAN)
    {
      Bool constrain = FALSE;
      Gnss::DesignMatrix A(normalEquationInfo, Vector(1));
      for(auto &trans : transmitter)
        for(auto &model : trans->biasModel)
          if(model.indexParameter)
          {
            Vector mean(model.temporal->parameterCount());
            for(const Time &time : times)
              mean += model.temporal->factors(time);

            const Double sigma = times.size()*0.0001; // 0.1 mm
            A.init(-1./sigma * mean.trans() * model.x); // constrain towards zero (0-x0)
            axpy(1./sigma, mean.trans(), A.column(model.indexParameter));
            A.accumulateNormals(normals, n, lPl, obsCount);
            constrain = TRUE;
          }
      if(constrain)
        logStatus<<"apply temporal zero mean of signal bias models of "<<name()<<" transmitters"<<Log::endl;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double GnssParametrizationTransmitter::updateParameter(const Gnss::NormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz, Bool printStatistics)
{
  try
  {
    // clock error
    // -----------
    Double      maxChangeClock = 0;
    std::string infoMaxChangeClock;
    for(auto &trans : transmitter)
      for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
        if(trans->indexParameterClock.at(idEpoch))
        {
          const Double dClock = x(normalEquationInfo.index(trans->indexParameterClock.at(idEpoch)), 0);
          trans->dclock(idEpoch) += dClock;

          if(std::fabs(dClock) > maxChangeClock)
          {
            maxChangeClock = std::fabs(dClock);
            infoMaxChangeClock = "  "+trans->name()+":  "+times.at(idEpoch).dateTimeStr()+" clockChange    = "+(1e3*dClock)%"%6.1f mm"s;
          }
        }

    // ============================================================

    // update positions & velocities
    // -----------------------------
    Double      maxChangePos = 0;
    std::string infoMaxChangePos;
    for(auto &trans : transmitter)
      if(trans->indexParameterSatellite && trans->indexParameterSatelliteArc)
      {
        // update transmitter xOrbit (estimated parameters vector)
        axpy(1., x.row(normalEquationInfo.index(trans->indexParameterSatellite), trans->countParameterSatellite),
             trans->xOrbit.row(0, trans->countParameterSatellite));
        axpy(1., x.row(normalEquationInfo.index(trans->indexParameterSatelliteArc), trans->countParameterSatelliteArc),
             trans->xOrbit.row(trans->countParameterSatellite, trans->countParameterSatelliteArc));

        Vector dpos = trans->variationalEquation.pos0;

        // gravity field parameters
        UInt countParameterGravityField = 0;
        if(gnss().gravityField && gnss().gravityField->normalEquationIndex())
        {
          countParameterGravityField = gnss().gravityField->parametrization()->parameterCount();
          matMult(1., trans->variationalEquation.PosDesign.column(0, countParameterGravityField),
                  x.row(normalEquationInfo.index(gnss().gravityField->normalEquationIndex()), countParameterGravityField),
                  trans->variationalEquation.pos0);
          matMult(1., trans->variationalEquation.VelDesign.column(0, countParameterGravityField),
                  x.row(normalEquationInfo.index(gnss().gravityField->normalEquationIndex()), countParameterGravityField),
                  trans->variationalEquation.vel0);
        }

        // variational equation parameters
        matMult(1., trans->variationalEquation.PosDesign.column(countParameterGravityField, trans->countParameterSatellite),
                x.row(normalEquationInfo.index(trans->indexParameterSatellite), trans->countParameterSatellite),
                trans->variationalEquation.pos0);
        matMult(1., trans->variationalEquation.VelDesign.column(countParameterGravityField, trans->countParameterSatellite),
                x.row(normalEquationInfo.index(trans->indexParameterSatellite), trans->countParameterSatellite),
                trans->variationalEquation.vel0);
        matMult(1., trans->variationalEquation.PosDesign.column(countParameterGravityField + trans->countParameterSatellite, trans->countParameterSatelliteArc),
                x.row(normalEquationInfo.index(trans->indexParameterSatelliteArc), trans->countParameterSatelliteArc),
                trans->variationalEquation.pos0);
        matMult(1., trans->variationalEquation.VelDesign.column(countParameterGravityField + trans->countParameterSatellite, trans->countParameterSatelliteArc),
                x.row(normalEquationInfo.index(trans->indexParameterSatelliteArc), trans->countParameterSatelliteArc),
                trans->variationalEquation.vel0);

        // max change of position
        dpos -= trans->variationalEquation.pos0;
        for(UInt i=0; i<trans->variationalEquation.times.size(); i++)
          if(norm(dpos.row(3*i,3)) > maxChangePos)
          {
            maxChangePos = norm(dpos.row(3*i,3));
            infoMaxChangePos = "  "+trans->name()+":  "+trans->variationalEquation.times.at(i).dateTimeStr()+" posChange      = "+(1e3*maxChangePos)%"%6.1f mm"s;
          }
      } // for(idTrans)

    // ============================================================

    // update signal biases
    // --------------------
    Double      maxChangeBias = 0;
    std::string infoMaxChangeBias;
    for(auto &trans : transmitter)
      if(trans->indexParameterBias)
      {
        const Vector dBias = trans->Bias * x.row(normalEquationInfo.index(trans->indexParameterBias), trans->Bias.columns());
        for(UInt idType=0; idType<dBias.size(); idType++)
          trans->signalBias.bias.at(idType) += dBias(idType);
        for(UInt idType=0; idType<dBias.size(); idType++)
          if(std::fabs(dBias(idType)) > maxChangeBias)
          {
            maxChangeBias = std::fabs(dBias(idType));
            infoMaxChangeBias = "  "+trans->name()+":  "+trans->signalBias.type.at(idType).str()+"              biasChange     = "+(1e3*dBias(idType))%"%6.1f mm"s;
          }
      } // for(idTrans)

    // ============================================================

    // update signal bias model
    // ------------------------
    Double      maxChangeBiasModel = 0;
    std::string infoMaxChangeBiasModel;
    for(auto &trans : transmitter)
      for(auto &model : trans->biasModel)
        if(model.indexParameter)
        {
          const Vector dx = x.row(normalEquationInfo.index(model.indexParameter), model.temporal->parameterCount());
          model.x += dx;
          for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
          {
            const Double dBias = inner(model.temporal->factors(times.at(idEpoch)), dx);
            model.bias(idEpoch) += dBias;
            if(std::fabs(dBias) > maxChangeBiasModel)
            {
              maxChangeBiasModel     = std::fabs(dBias);
              infoMaxChangeBiasModel = "  "+trans->name()+":  "+times.at(idEpoch).dateTimeStr()+" biasModel      = "+(1e3*dBias)%"%6.1f mm"s;
            }
          }
        } // for(trans, model)

    // ============================================================

    // update clock model
    // ------------------
    Double      maxChangeClockModel = 0;
    std::string infoMaxChangeClockModel;
    for(auto &trans : transmitter)
      if(trans->indexParameterClockModel)
      {
        const Vector dx = x.row(normalEquationInfo.index(trans->indexParameterClockModel), clockModel->parameterCount());
        for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
        {
          const Double dClock = inner(clockModel->factors(times.at(idEpoch)), dx);
          trans->clock0.at(idEpoch) += dClock;
          if(std::fabs(dClock) > maxChangeClockModel)
          {
            maxChangeClockModel     = std::fabs(dClock);
            infoMaxChangeClockModel = "  "+trans->name()+":  "+times.at(idEpoch).dateTimeStr()+" clockModel     = "+(1e3*dClock)%"%6.1f mm"s;
          }
        }
      } // for(idTrans)

    // ============================================================

    Gnss::checkMaxChange(maxChangePos,         infoMaxChangePos,         printStatistics, normalEquationInfo.comm);
    Gnss::checkMaxChange(maxChangeClock,       infoMaxChangeClock,       printStatistics, normalEquationInfo.comm);
    Gnss::checkMaxChange(maxChangeClockModel,  infoMaxChangeClockModel,  printStatistics, normalEquationInfo.comm);
    Gnss::checkMaxChange(maxChangeBias,        infoMaxChangeBias,        printStatistics, normalEquationInfo.comm);
    Gnss::checkMaxChange(maxChangeBiasModel,   infoMaxChangeBiasModel,   printStatistics, normalEquationInfo.comm);

    // ============================================================

    // VCE for clock noise
    // -------------------
    if(Wz.size() && estimateClockSigma && (estimateClockError == EstimateClockError::NOISEMODEL))
    {
      for(auto &trans : transmitter)
        if(trans->useable())
        {
          Double ePe        = 0;
          Double redundancy = 0;
          Double obsCount   = 0;
          for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
            if(trans->indexParameterClock.at(idEpoch))
            {
              // how many epochs back in time?
              UInt count = 0;
              while(count<trans->clockNoise.rows() && (count<=idEpoch) && trans->indexParameterClock.at(idEpoch-count))
                count++;

              Double We = 0;
              Matrix Awz(1, Wz.columns());
              for(UInt i=0; i<count; i++)
              {
                We += trans->clockNoise(count-1, count-1-i) * trans->dclock(idEpoch-i);
                axpy(trans->clockNoise(count-1, count-1-i), Wz.row(normalEquationInfo.index(trans->indexParameterClock.at(idEpoch-i))), Awz);
              }

              ePe        += We*We;
              redundancy += 1 - quadsum(Awz);
              obsCount++;
            }

          if(obsCount)
          {
            const Double factor = sqrt(ePe/redundancy);
            trans->clockSigma *= factor;
            trans->clockNoise *= 1./factor;
            if(printStatistics)
              logInfo<<"  "<<trans->name()<<": sigma clock noise = "<<trans->clockSigma%"%4.2f redundancy = "s<<(redundancy/obsCount)%"%4.2f count = "s<<obsCount%"%5i"s<<Log::endl;
          }
        }
    }

    // ============================================================

    return std::max(maxChangeClock, std::max(maxChangeClockModel, std::max(maxChangePos, std::max(maxChangeBias, maxChangeBiasModel))));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationTransmitter::writeResults(const Gnss::NormalEquationInfo &normalEquationInfo, const std::string &suffix)
{
  try
  {
    if(!Parallel::isMaster(normalEquationInfo.comm))
      return;

    VariableList fileNameVariableList;
    addVariable("prn", fileNameVariableList);

    // write used transmitter list
    // ---------------------------
    if(!fileNameOutTransmitterList.empty() && Parallel::isMaster(normalEquationInfo.comm))
    {
      logStatus<<"write used transmitter list to file <"<<fileNameOutTransmitterList(fileNameVariableList).appendBaseName(suffix)<<">"<<Log::endl;
      std::vector<std::string> usedTransmitterList;
      for(const auto &trans : transmitter)
        if(trans->useable())
          usedTransmitterList.push_back(trans->name());
      writeFileStringList(fileNameOutTransmitterList(fileNameVariableList).appendBaseName(suffix), usedTransmitterList);
    }

    // variational equations and orbits
    // --------------------------------
    if(estimateDynamicOrbit && (!fileNameTemplateOutVariational.empty() || !fileNameTemplateOutOrbit.empty()))
    {
      fileNameVariableList["prn"]->setValue("***");
      if(!fileNameTemplateOutOrbit.empty())
        logStatus<<"write transmitter orbits to files <"<<fileNameTemplateOutOrbit(fileNameVariableList).appendBaseName(suffix)<<">"<<Log::endl;
      if(!fileNameTemplateOutVariational.empty())
        logStatus<<"write variational equations to files <"<<fileNameTemplateOutVariational(fileNameVariableList).appendBaseName(suffix)<<">"<<Log::endl;

      UInt countParameterGravityField = 0;
      if(gnss().gravityField)
        countParameterGravityField = gnss().gravityField->parametrization()->parameterCount();

      for(auto &trans : transmitter)
        if(trans->useable())
        {
          fileNameVariableList["prn"]->setValue(trans->name());

          VariationalEquationArc variationalArc;
          variationalArc.times = trans->variationalEquation.times;
          variationalArc.pos0  = trans->variationalEquation.pos0;
          variationalArc.vel0  = trans->variationalEquation.vel0;

          if(!fileNameTemplateOutVariational.empty())
          {
            // copy state transition matrix and satellite orientation from input variational file
            variationalArc.PosState = trans->variationalEquation.PosDesign.column(countParameterGravityField + trans->countParameterSatellite, trans->countParameterSatelliteArc);
            variationalArc.VelState = trans->variationalEquation.VelDesign.column(countParameterGravityField + trans->countParameterSatellite, trans->countParameterSatelliteArc);
            variationalArc.rotSat   = trans->variationalEquation.rotSat;

            // (updated) earth rotation
            variationalArc.rotEarth.resize(variationalArc.times.size());
            for(UInt idEpoch=0; idEpoch<variationalArc.times.size(); idEpoch++)
              variationalArc.rotEarth.at(idEpoch) = gnss().earthRotation->rotaryMatrix(variationalArc.times.at(idEpoch));

            writeFileVariationalEquation(fileNameTemplateOutVariational(fileNameVariableList).appendBaseName(suffix), trans->satelliteModel, variationalArc);
          }

          if(!fileNameTemplateOutOrbit.empty())
            InstrumentFile::write(fileNameTemplateOutOrbit(fileNameVariableList).appendBaseName(suffix), variationalArc.orbitArc());
        }
    }

    // Estimated parameters
    // --------------------
    if(estimateDynamicOrbit && !fileNameTemplateOutParameter.empty())
    {
      fileNameVariableList["prn"]->setValue("***");
      logStatus<<"write estimated transmitter parameters to files <"<<fileNameTemplateOutParameter(fileNameVariableList).appendBaseName(suffix)<<">"<<Log::endl;
      for(auto &trans : transmitter)
        if(trans->useable())
        {
          fileNameVariableList["prn"]->setValue(trans->name());
          writeFileMatrix(fileNameTemplateOutParameter(fileNameVariableList).appendBaseName(suffix), trans->xOrbit);
        }
    }

    // Signal bias
    // -----------
    if(!fileNameTemplateOutSignal.empty())
    {
      fileNameVariableList["prn"]->setValue("***");
      logStatus<<"write transmitter signal bias to files <"<<fileNameTemplateOutSignal(fileNameVariableList).appendBaseName(suffix)<<">"<<Log::endl;
      for(auto &trans : transmitter)
        if(trans->useable())
        {
          GnssSignalBias signalBias = trans->signalBias;
          for(UInt idType=0; idType<signalBias.type.size(); idType++)
            if(signalBias.type.at(idType) == GnssType::PHASE)
              signalBias.bias.at(idType) = std::remainder(signalBias.bias.at(idType), LIGHT_VELOCITY/signalBias.type.at(idType).frequency());
          fileNameVariableList["prn"]->setValue(trans->name());
          writeFileGnssSignalBias(fileNameTemplateOutSignal(fileNameVariableList).appendBaseName(suffix), signalBias);
        }
    }

    // Signal bias model
    // -----------------
    for(auto &trans : transmitter)
      if(trans->useable())
        for(auto &model : trans->biasModel)
          if(!model.fileNameOut.empty())
          {
            fileNameVariableList["prn"]->setValue(trans->name());
            MiscValueArc arc;
            for(UInt idEpoch : normalEquationInfo.idEpochs)
              if(trans->useable(idEpoch))
              {
                const Double bias = model.bias(idEpoch);
                if(bias)
                {
                  MiscValueEpoch epoch;
                  epoch.time  = times.at(idEpoch);
                  epoch.value = bias;
                  arc.push_back(epoch);
                }
              }

            if(arc.size())
              InstrumentFile::write(model.fileNameOut(fileNameVariableList).appendBaseName(suffix), arc);
          }

    // Estimated clocks
    // ----------------
    if((estimateClockError != EstimateClockError::NONE) && !fileNameOutClock.empty())
    {
      fileNameVariableList["prn"]->setValue("***");
      logStatus<<"write transmitter clocks to files <"<<fileNameOutClock(fileNameVariableList).appendBaseName(suffix)<<">"<<Log::endl;
      for(auto &trans : transmitter)
        if(trans->useable())
        {
          fileNameVariableList["prn"]->setValue(trans->name());

          MiscValueArc arc;
          for(UInt idEpoch : normalEquationInfo.idEpochs)
            if(trans->useable(idEpoch))
            {
              const Double clk = trans->clockError(idEpoch, times.at(idEpoch));
              if(clk)
              {
                MiscValueEpoch epoch;
                epoch.time  = times.at(idEpoch);
                epoch.value = clk;
                arc.push_back(epoch);
              }
            }

          if(arc.size())
            InstrumentFile::write(fileNameOutClock(fileNameVariableList).appendBaseName(suffix), arc);
        }
    } // if(estimateClockError)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

Bool GnssParametrizationTransmitter::Transmitter::isDesignMatrixTransmitter(const Gnss::NormalEquationInfo &/*normalEquationInfo*/, UInt /*idRecv*/, UInt idEpoch) const
{
  try
  {
    if(indexParameterClock.at(idEpoch) ||
       indexParameterSatellite         ||
       indexParameterSatelliteArc      ||
       indexParameterBias              ||
       indexParameterClockModel)
      return TRUE;

    for(auto &model : biasModel)
      if(model.indexParameter)
        return TRUE;

    return FALSE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationTransmitter::Transmitter::designMatrixTransmitter(const Gnss::NormalEquationInfo &/*normalEquationInfo*/, const Gnss::ObservationEquation &eqn, Gnss::DesignMatrix &A) const
{
  try
  {
    // clock error
    // -----------
    if(indexParameterClock.at(eqn.idEpoch))
      copy(eqn.A.column(Gnss::ObservationEquation::idxClockTrans,1), A.column(indexParameterClock.at(eqn.idEpoch)));

    // ===========================================================

    // variational equations
    // ---------------------
    if(indexParameterSatellite && indexParameterSatelliteArc)
    {
      // gravity field
      UInt countParameterGravityField = 0;
      if(gnss().gravityField && gnss().gravityField->normalEquationIndex())
      {
        countParameterGravityField = gnss().gravityField->parametrization()->parameterCount();
        matMult(1., eqn.A.column(Gnss::ObservationEquation::idxPosTrans,3),
                base->polynomial.interpolate({eqn.timeTrans}, variationalEquation.times, variationalEquation.PosDesign.column(0, countParameterGravityField), 3),
                A.column(gnss().gravityField->normalEquationIndex()));
      } // if(base->countParameterGravityField)

      // variational equations
      matMult(1., eqn.A.column(Gnss::ObservationEquation::idxPosTrans,3),
              base->polynomial.interpolate({eqn.timeTrans}, variationalEquation.times, variationalEquation.PosDesign.column(countParameterGravityField, countParameterSatellite), 3),
              A.column(indexParameterSatellite));
      matMult(1., eqn.A.column(Gnss::ObservationEquation::idxPosTrans,3),
              base->polynomial.interpolate({eqn.timeTrans}, variationalEquation.times, variationalEquation.PosDesign.column(countParameterGravityField + countParameterSatellite, countParameterSatelliteArc), 3),
              A.column(indexParameterSatelliteArc));
    }

    // ===========================================================

    const Time timeTrans = std::max(eqn.timeTrans, base->times.at(0));

    // signal biases
    // -------------
    if(indexParameterBias)
    {
      MatrixSlice Design(A.column(indexParameterBias));
      for(UInt idType=0; idType<eqn.typesTransmitted.size(); idType++)
      {
        const UInt idx = GnssType::index(signalBias.type, eqn.typesTransmitted.at(idType));
        if(idx != NULLINDEX)
          matMult(1., eqn.A.column(Gnss::ObservationEquation::idxUnit + eqn.types.size() + idType), Bias.row(idx), Design);
      }
    }

    for(auto &model : biasModel)
      if(model.indexParameter)
      {
        MatrixSlice Design(A.column(model.indexParameter));
        for(UInt idType=0; idType<eqn.typesTransmitted.size(); idType++)
          if(eqn.typesTransmitted.at(idType) == model.type)
            matMult(1., eqn.A.column(Gnss::ObservationEquation::idxUnit + eqn.types.size() + idType), model.temporal->factors(timeTrans).trans(), Design);
      }

    // ===========================================================

    // clock model
    // -----------
    if(indexParameterClockModel)
    {
      matMult(1., eqn.A.column(Gnss::ObservationEquation::idxClockTrans,1), base->clockModel->factors(timeTrans).trans(),
              A.column(indexParameterClockModel));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool GnssParametrizationTransmitter::Transmitter::supportsIntegerAmbiguities(const Gnss::NormalEquationInfo &/*normalEquationInfo*/) const
{
  return base->integerAmbiguities;
}

/***********************************************/

Bool GnssParametrizationTransmitter::Transmitter::isCodeBiasEstimated(const Gnss::NormalEquationInfo &normalEquationInfo) const
{
  try
  {
    return base->estimateCodeBias && (normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_TRANSMITTER_SIGNALBIAS);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool GnssParametrizationTransmitter::Transmitter::isPhaseBiasEstimated(const Gnss::NormalEquationInfo &normalEquationInfo) const
{
  try
  {
    return base->estimatePhaseBias && (normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_TRANSMITTER_SIGNALBIAS);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string GnssParametrizationTransmitter::Transmitter::name() const
{
  return PRN().str().substr(3,3);
}

/***********************************************/

Bool GnssParametrizationTransmitter::Transmitter::useable(UInt idEpoch) const
{
  if(!useAtAll)
    return FALSE;
  if(idEpoch == NULLINDEX)
    return useAtAll;
  return use.at(idEpoch);
}

/***********************************************/

void GnssParametrizationTransmitter::Transmitter::disable(UInt idEpoch)
{
  if(idEpoch == NULLINDEX)
    useAtAll = FALSE;
  else
    use.at(idEpoch) = FALSE;
}

/***********************************************/

Vector3d GnssParametrizationTransmitter::Transmitter::positionCoM(UInt /*idEpoch*/, const Time &time) const
{
  try
  {
    Vector p = base->polynomial.interpolate({time}, variationalEquation.times, variationalEquation.pos0, 3);
    return Vector3d(p(0), p(1), p(2));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector3d GnssParametrizationTransmitter::Transmitter::position(UInt idEpoch, const Time &time) const
{
  try
  {
    const UInt idAntenna = transmitterInfo.findAntenna(time);
    if(idAntenna >= transmitterInfo.antenna.size())
      throw(Exception(name() + ": no antenna found at " + time.dateTimeStr()));
    const Vector3d offset = transmitterInfo.antenna.at(idAntenna).position - transmitterInfo.referencePoint(time); // anntenna center relative to center of mass (CoM) in body system
    return positionCoM(idEpoch, time) + celestial2antennaFrame(idEpoch, time).inverseTransform(transmitterInfo.antenna.at(idAntenna).local2antennaFrame.transform(offset));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector3d GnssParametrizationTransmitter::Transmitter::velocity(UInt /*idEpoch*/, const Time &time) const
{
  try
  {
    Vector v = base->polynomial.interpolate({time}, variationalEquation.times, variationalEquation.vel0, 3);
    return Vector3d(v(0), v(1), v(2));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Transform3d GnssParametrizationTransmitter::Transmitter::celestial2antennaFrame(UInt idEpoch, const Time &/*time*/) const
{
  return crf2arf.at(idEpoch);
}

/***********************************************/

Double GnssParametrizationTransmitter::Transmitter::clockError(UInt idEpoch, const Time &time) const
{
  try
  {
    // return value without interpolation if time is close to times.at(idEpoch)
    if(std::fabs((base->times.at(idEpoch)-time).seconds()) < 0.1)
      return (clock0(idEpoch)+dclock(idEpoch))/LIGHT_VELOCITY;

    // switch to previous interval for interpolation if time < times.at(idEpoch)
    if((base->times.at(idEpoch)>time && idEpoch!=0) || idEpoch==clock0.size()-1)
      idEpoch--;

    // linear interpolation
    const Double t = (time-base->times.at(idEpoch)).seconds()/(base->times.at(idEpoch+1)-base->times.at(idEpoch)).seconds();
    return ((1-t)*(clock0(idEpoch)+dclock(idEpoch)) + t*(clock0(idEpoch+1)+dclock(idEpoch+1)))/LIGHT_VELOCITY;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::vector<GnssType> GnssParametrizationTransmitter::Transmitter::definedTypes(UInt idEpoch) const
{
  try
  {
    UInt idRecv = transmitterInfo.findReceiver(base->times.at(idEpoch));
    if(idRecv != NULLINDEX && transmitterInfo.receiver.at(idRecv).receiverDef)
      return transmitterInfo.receiver.at(idRecv).receiverDef->types;

    return std::vector<GnssType>();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector GnssParametrizationTransmitter::Transmitter::antennaVariations(UInt idEpoch, Angle azimut, Angle elevation, const std::vector<GnssType> &type) const
{
  try
  {
    Vector corr(type.size());
    corr += transmitterInfo.antennaVariations(base->times.at(idEpoch), azimut, elevation, type, base->noPatternFoundAction);
    corr += signalBias.compute(type); // Code/Phase biases

    for(auto &model : biasModel)
      for(UInt idType=0; idType<type.size(); idType++)
        if(type.at(idType) == model.type)
          corr(idType) += model.bias(idEpoch);

    for(auto &bias : timeVariableBiases)
      for(UInt idType=0; idType<type.size(); idType++)
        if(type.at(idType) == bias.type)
          corr(idType) += bias.biases(idEpoch);

    return corr;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
