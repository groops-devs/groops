/***********************************************/
/**
* @file gnssParametrizationReceiverLeo.cpp
*
* @brief GNSS for Low Earth Orbiter (LEO).
*
* @author Torsten Mayer-Guerr
* @author Norbert Zehentner
* @author Sebastian Strasser
* @date 2010-08-03
*
*/
/***********************************************/

#include "base/import.h"
#include "parallel/matrixDistributed.h"
#include "config/config.h"
#include "inputOutput/logging.h"
#include "parser/dataVariables.h"
#include "parser/expressionParser.h"
#include "files/fileInstrument.h"
#include "files/fileGnssStationInfo.h"
#include "files/fileGnssSignalBias.h"
#include "files/fileParameterName.h"
#include "gnss/gnss.h"
#include "gnss/gnssObservation.h"
#include "gnss/gnssDesignMatrix.h"
#include "gnss/gnssTransmitter.h"
#include "gnss/gnssParametrizationReceiverLeo.h"

/***********************************************/

GnssParametrizationReceiverLeo::GnssParametrizationReceiverLeo(Config &config)
{
  try
  {
    noPatternFoundAction = GnssAntennaDefinition::NoPatternFoundAction::IGNORE_OBSERVATION;
    FileName fileNameStationInfo, fileNameAntennaDef, fileNameReceiverDef, fileNameAccuracyDef;

    readConfig(config, "outputfileOrbit",              fileNameOutOrbit,         Config::OPTIONAL, "output/orbit_{loopTime:%D}.dat",      "estimated kinematic orbit");
    readConfig(config, "outputfileCovariancePodEpoch", fileNameOutCovariance,    Config::OPTIONAL, "output/covariance_{loopTime:%D}.dat", "3x3 epoch covariances");
    readConfig(config, "outputfileClock",              fileNameOutClock,         Config::OPTIONAL, "output/clock_{loopTime:%D}.dat",      "estimtated clock errors");
    readConfig(config, "outputfileResiduals",          fileNameResiduals,        Config::OPTIONAL, "output/residuals_{loopTime:%D}.dat",  "");
    readConfig(config, "outputfileSignalBias",         fileNameOutSignal,        Config::OPTIONAL, "output/signaBias_{loopTime:%D}.txt",  "estimated signal biases");
    readConfig(config, "inputfileStationInfo",         fileNameStationInfo,      Config::MUSTSET,  "", "satellite metadata (antenna, receiver, ...)");
    readConfig(config, "inputfileAntennaDefinition",   fileNameAntennaDef,       Config::MUSTSET,  "", "antenna center offsets and variations");
    readConfig(config, "inputfileReceiverDefinition",  fileNameReceiverDef,      Config::OPTIONAL, "", "observed signal types");
    readConfig(config, "inputfileAccuracyDefinition",  fileNameAccuracyDef,      Config::MUSTSET,  "", "elevation and azimut dependent accuracy");
    readConfig(config, "inputfileSignalBias",          fileNameSignalBias,       Config::OPTIONAL, "", "a priori signal biases");
    readConfig(config, "inputfileObservations",        fileNameObs,              Config::MUSTSET,  "gnssReceiver_{loopTime:%D}.dat", "");
    readConfig(config, "inputfileOrbit",               fileNameOrbit,            Config::MUSTSET,  "",    "approximate positions");
    readConfig(config, "inputfileStarCamera",          fileNameStarCamera,       Config::MUSTSET,  "",    "satellite attitude");
    readConfig(config, "estimateCodeBias",             estimateCodeBias,         Config::DEFAULT,  "0",   "");
    readConfig(config, "estimatePhaseBias",            estimatePhaseBias,        Config::DEFAULT,  "1",   "");
    readConfig(config, "supportsIntegerAmbiguities",   integerAmbiguities,       Config::DEFAULT,  "1",   "receiver tracks full cycle integer ambiguities");
    readConfig(config, "wavelengthFactor",             wavelengthFactor_,        Config::DEFAULT,  "1.",  "factor to account for half-wavelength observations (collected by codeless squaring techniques)");
    readConfig(config, "sigmaFactorPhase",             exprSigmaPhase,           Config::OPTIONAL, "",    "PHASE: factor = f(FREQ, ELE, SNR, ROTI, dTEc, IONOINDEX)");
    readConfig(config, "sigmaFactorCode",              exprSigmaCode,            Config::OPTIONAL, "",    "CODE: factor = f(FREQ, ELE, SNR, ROTI, dTEc, IONOINDEX)");
    readConfig(config, "antennaCenter",                antennaCenterVariations,  Config::OPTIONAL, "",    "estimate antenna center variations");
    if(readConfigSequence(config, "preprocessing", Config::MUSTSET, "", ""))
    {
      readConfig(config, "elevationCutOff",       cutOff,                Config::DEFAULT,  "0",   "[degree] ignore observations below cutoff");
      readConfig(config, "elevationTrackMinimum", elevationTrackMinimum, Config::DEFAULT,  "0",   "[degree] ignore tracks that never exceed minimum elevation");
      readConfig(config, "useType",               useType,               Config::OPTIONAL, "" ,   "only use observations that match any of these patterns");
      readConfig(config, "ignoreType",            ignoreType,            Config::OPTIONAL, "",    "ignore observations that match any of these patterns");
      std::string choice;
      if(readConfigChoice(config, "noAntennaPatternFound", choice, Config::MUSTSET, "ignoreObservation", "what should happen if no antenna pattern is found for an observation"))
      {
        if(readConfigChoiceElement(config, "ignoreObservation",   choice, "ignore observation if no matching pattern is found"))
          noPatternFoundAction = GnssAntennaDefinition::NoPatternFoundAction::IGNORE_OBSERVATION;
        if(readConfigChoiceElement(config, "useNearestFrequency", choice, "use pattern of nearest frequency if not matching pattern is found"))
          noPatternFoundAction = GnssAntennaDefinition::NoPatternFoundAction::USE_NEAREST_FREQUENCY;
        if(readConfigChoiceElement(config, "throwException",      choice, "throw exception if no matching pattern is found"))
          noPatternFoundAction = GnssAntennaDefinition::NoPatternFoundAction::THROW_EXCEPTION;
        endChoice(config);
      }
      readConfig(config, "huber",                 huber,                 Config::DEFAULT,  "2.5", "residuals > huber*sigma0 are downweighted");
      readConfig(config, "huberPower",            huberPower,            Config::DEFAULT,  "1.5", "residuals > huber: sigma=(e/huber)^huberPower*sigma0");
      readConfig(config, "codeMaxPositionDiff",   codeMaxPosDiff,        Config::DEFAULT,  "100", "[m] max. allowed position error by PPP code only clock error estimation");
      readConfig(config, "minObsCountPerTrack",   minObsCountPerTrack,   Config::DEFAULT,  "20",  "tracks with less number of epochs with observations are dropped");
      readConfig(config, "windowROTI",            windowROTI,            Config::DEFAULT,  "0",   "[seconds] moving time window to compute Rate Of Tec Index (ROTI) for each observation");
      endSequence(config);
    } // readConfigSequence(preprocessing)
    if(isCreateSchema(config)) return;

    // station info and antenna definition
    // -----------------------------------
    std::vector<GnssAntennaDefinitionPtr> antennaDefList;
    if(!fileNameAntennaDef.empty())
      readFileGnssAntennaDefinition(fileNameAntennaDef, antennaDefList);

    std::vector<GnssReceiverDefinitionPtr> receiverDefList;
    if(!fileNameReceiverDef.empty())
      readFileGnssReceiverDefinition(fileNameReceiverDef, receiverDefList);

    std::vector<GnssAntennaDefinitionPtr> accuracyDefList;
    if(!fileNameAccuracyDef.empty())
      readFileGnssAntennaDefinition(fileNameAccuracyDef, accuracyDefList);

    readFileGnssStationInfo(fileNameStationInfo, stationInfo);
    stationInfo.fillAntennaPattern (antennaDefList);
    stationInfo.fillReceiverDefinition(receiverDefList);
    stationInfo.fillAntennaAccuracy(accuracyDefList);

    // create variables for sigma computation
    // --------------------------------------
    varList = config.getVarList();
    addVariable("SNR",       varList);
    addVariable("FREQ",      varList);
    addVariable("ROTI",      varList);
//     addVariable("ELE",       varList);
//     addVariable("dTEC",      varList);
//     addVariable("IONOINDEX", varList);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::vector<Gnss::ParametrizationPtr> GnssParametrizationReceiverLeo::parametrizations()
{
  std::vector<Gnss::ParametrizationPtr> p({shared_from_this()});
  p.insert(p.begin(), antennaCenterVariations.begin(), antennaCenterVariations.end());
  return p;
}

/***********************************************/

void GnssParametrizationReceiverLeo::initIntervalReceiver(Gnss::AnalysisType analysisType, const std::vector<Time> &times, const Time &timeMargin, Parallel::CommunicatorPtr comm)
{
  try
  {
    // delete old observations
    free();
    pos0.clear();
    displacement.clear();
    dpos.clear();
    cov.clear();
    vel.clear();
    crf2sat.clear();
    indexParameterEpoch.clear();
    indexParameterEpoch.resize(times.size(), Gnss::ParameterIndex());
    signalBias = GnssSignalBias();

    // ===========================================================

    // init
    // ----
    if(Parallel::isMaster(comm))
    {
      initInterval(times);
      pos0.resize(times.size());
      displacement.resize(times.size());
      dpos.resize(times.size());
      cov.resize(times.size());
      vel.resize(times.size());
      crf2sat.resize(times.size());
    }

    timeIntervalStart = times.front();
    timeIntervalEnd   = times.back()+medianSampling(times);
    VariableList fileNameVariableList;
    addTimeVariables(fileNameVariableList);
    evaluateTimeVariables(0, timeIntervalStart, timeIntervalEnd, fileNameVariableList);

    // test completeness of antennas
    for(const auto &antenna : stationInfo.antenna)
      if(antenna.timeEnd >= times.front() && antenna.timeStart < times.back() && (!antenna.antennaDef || !antenna.accuracyDef))
      {
        logWarning<<stationInfo.markerName<<"."<<stationInfo.markerNumber<<": no "<<(!antenna.antennaDef ? "antenna" : "accuracy")
                  <<" definition found for "<<antenna.str()<<", disabling receiver."<<Log::endl;
        disable();
        return;
      }

    // read approximate positions
    // --------------------------
    if(useable())
    {
      try
      {
        logStatus<<"read approx. orbit <"<<fileNameOrbit(fileNameVariableList)<<">"<<Log::endl;
        InstrumentFile fileOrbit(fileNameOrbit(fileNameVariableList));
        InstrumentFile fileStarCamera(fileNameStarCamera(fileNameVariableList));
        InstrumentFile::checkArcCount({fileOrbit, fileStarCamera});

        UInt idEpoch = 0;
        for(UInt arcNo=0; arcNo<fileOrbit.arcCount(); arcNo++)
        {
          OrbitArc      orbit      = fileOrbit.readArc(arcNo);
          StarCameraArc starCamera = fileStarCamera.readArc(arcNo);
          Arc::checkSynchronized({orbit, starCamera});
          for(UInt i=0; i<orbit.size(); i++)
          {
            while((idEpoch < times.size()) && (times.at(idEpoch) < orbit.at(i).time-timeMargin))
              disable(idEpoch++);
            if(idEpoch >= times.size())
              break;
            if(times.at(idEpoch) > orbit.at(i).time+timeMargin)
              continue;

            // this->times.at(idEpoch) = orbit.at(i).time; // overwritten by readObservations()
            this->pos0.at(idEpoch)  = orbit.at(i).position;
            this->vel.at(idEpoch)   = orbit.at(i).velocity;

            if(pos0.at(idEpoch).r() == 0.)
            {
              logWarning<<"pos.at(idEpoch).r() == 0 at "<<orbit.at(i).time.dateTimeStr()<<Log::endl;
              disable(idEpoch);
              continue;
            }

            // calculate antenna center in CRF
            // -------------------------------
            const UInt idAnt = stationInfo.findAntenna(timeCorrected(idEpoch));
            if(idAnt == NULLINDEX)
              throw(Exception(stationInfo.markerName+"."+stationInfo.markerNumber+": no antenna definition found at "+timeCorrected(idEpoch).dateTimeStr()));

            // orientation of antenna
            this->crf2sat.at(idEpoch) = inverse(starCamera.at(i).rotary);
            // correct for antenna position
            this->displacement.at(idEpoch) = starCamera.at(i).rotary.rotate(stationInfo.antenna.at(idAnt).position - stationInfo.referencePoint(timeCorrected(idEpoch)));

            idEpoch++;
          }
          if(idEpoch >= times.size())
            break;
        }
      }
      catch(std::exception &e)
      {
        logWarning<<name()<<": Initialization of satellite failed, disabling. "<<e.what()<< Log::endl;
        disable();
        return;
      }
    }

    // ===========================================================

    // code & phase biases
    // -------------------
    if(!fileNameSignalBias.empty())
      readFileGnssSignalBias(fileNameSignalBias(fileNameVariableList), signalBias);

    // ===========================================================

    // read observations
    // -----------------
    if(!fileNameObs.empty() && useable())
    {
      try
      {
        logStatus<<"read GNSS receiver data <"<<fileNameObs(fileNameVariableList)<<">"<<Log::endl;
        InstrumentFile fileReceiver(fileNameObs(fileNameVariableList));
        readObservations(fileReceiver, useType, ignoreType, timeMargin);
      }
      catch(std::exception &e)
      {
        logWarning<<name()<<": Initialization of satellite failed, disabling. "<<e.what()<< Log::endl;
        disable();
        return;
      }
    }

    // ===========================================================

    // init observations
    // -----------------
    if(!fileNameObs.empty() && useable())
    {
      try
      {
        logStatus<<"init observations"<<Log::endl;
        initObservation(analysisType);
        deleteUndefinedObservations();
        if(analysisType & Gnss::ANALYSIS_PHASE)
          createTracks(minObsCountPerTrack, {GnssType::L5_G});
        deleteObservationsOfInestimableEpochs(analysisType);
        estimateInitialClockErrorFromCodeObservations(codeMaxPosDiff);

        ObservationEquationList eqn(Receiver::gnss(), *this, analysisType & (Gnss::ANALYSIS_CODE | Gnss::ANALYSIS_PHASE));
        disableEpochsWithGrossCodeObservationOutliers(eqn, codeMaxPosDiff);
        if(analysisType & Gnss::ANALYSIS_PHASE)
        {
          removeTracksWithInsufficientObservationCount(eqn, minObsCountPerTrack, analysisType);
          cycleSlipsDetection(eqn, 5/*lambda*/, 0/*tecWindowSize*/, 3.5/*tecSigmaFactor*/);
          // removeLowElevationTracks(eqn, elevationTrackMinimum);
          removeTracksWithInsufficientObservationCount(eqn, minObsCountPerTrack, analysisType);
          trackOutlierDetection(eqn, {GnssType::L5_G});
          cycleSlipsRepairAtSameFrequency(eqn);
          if(windowROTI > 0)
            addRotiPseudoObservations(eqn, windowROTI);
        }

        // apply factors for accuracies from expressions
        if(exprSigmaPhase || exprSigmaCode)
        {
          for(UInt idTrans=0; idTrans<Receiver::gnss().transmitter.size(); idTrans++)
            for(UInt idEpoch=0; idEpoch<idEpochSize(); idEpoch++)
            {
              Gnss::Observation *obs = observation(idTrans, idEpoch);
              if(obs)
              {
                varList["ROTI"]->setValue(0);
                const UInt idx = obs->index(GnssType::ROTI);
                if(idx != NULLINDEX)
                  varList["ROTI"]->setValue(obs->at(idx).observation);

                for(UInt idType=0; idType<obs->size(); idType++)
                  if((obs->at(idType).type == GnssType::RANGE) || (obs->at(idType).type == GnssType::PHASE))
                  {
                    varList["FREQ"]->setValue( obs->at(idType).type.frequency() );
                    varList["SNR"]->setValue(0);
                    const UInt idx = obs->index(GnssType::SNR + (obs->at(idType).type & GnssType::FREQUENCY));
                    if(idx != NULLINDEX)
                      varList["SNR"]->setValue(obs->at(idx).observation);

                    if(exprSigmaPhase && (obs->at(idType).type == GnssType::PHASE))
                    {
                      const Double factor = exprSigmaPhase->evaluate(varList);
                      obs->at(idType).sigma0 *= factor;
                      obs->at(idType).sigma  *= factor;
                    }

                    if(exprSigmaCode && (obs->at(idType).type == GnssType::RANGE))
                    {
                      const Double factor = exprSigmaCode->evaluate(varList);
                      obs->at(idType).sigma0 *= factor;
                      obs->at(idType).sigma  *= factor;
                    }
                  } // for(idType)
              } // if(obs)
            } // for(idTrans, idEpoch)
        } // if(exprSigmaPhase || exprSigmaCode)
      }
      catch(std::exception &e)
      {
        logWarning<<name()<<" disabled: "<<e.what()<<Log::endl;
        disable();
      }

      // test epochs
      UInt countEpoch = 0;
      for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
        if(useable(idEpoch))
          countEpoch++;
      if(countEpoch == 0)
      {
        logWarning<<name()<<" disabled: no valid epochs"<<Log::endl;
        disable();
      }
    } // init observations
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool GnssParametrizationReceiverLeo::initIntervalDisabling(Gnss::AnalysisType analysisType, const std::vector<Time> &times, const Time &/*timeMargin*/, Parallel::CommunicatorPtr comm)
{
  try
  {
    UInt disabled = FALSE;
    if(useable())
    {
      Gnss::Receiver::ObservationEquationList eqn;
      disabled = deleteObservationsOfInestimableEpochs(analysisType) || disabled;
      disabled = removeTracksWithInsufficientObservationCount(eqn, minObsCountPerTrack, analysisType) || disabled;

      // test epochs
      UInt countEpoch = 0;
      for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
        if(useable(idEpoch))
          countEpoch++;

      if(countEpoch == 0)
      {
        logWarning<<name()<<" disabled: no valid epochs"<<Log::endl;
        disable();
        disabled = TRUE;
      }
    }

    Parallel::reduceSum(disabled, 0, comm);
    Parallel::broadCast(disabled, 0, comm);
    return disabled;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationReceiverLeo::initIntervalLate(Gnss::AnalysisType /*analysisType*/, const std::vector<Time> &/*times*/, const Time &/*timeMargin*/, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    // get observation types
    std::vector<GnssType> types;
    for(auto &typesTrans : Receiver::gnss().typesRecvTrans.at(idRecv()))
      for(GnssType type : typesTrans)
        if((type == GnssType::PHASE) || (type == GnssType::RANGE))
          if(GnssType::index(types, type) == NULLINDEX)
            types.push_back(type);
    std::sort(types.begin(), types.end());

    // code & phase biases
    // -------------------
    if(types.size()) // consider simulation case
    {
      signalBias.bias = signalBias.compute(types); // apriori signal bias
      signalBias.type = types;
    }

    // antenna center variations
    // -------------------------
    for(auto &acv : antennaCenterVariations)
      acv->addAntennas(idRecv(), TRUE/*deviceIsReceiver*/, stationInfo, types);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void GnssParametrizationReceiverLeo::initParameter(Gnss::NormalEquationInfo &normalEquationInfo)
{
  try
  {
    // reset parameter indices
    // -----------------------
    std::fill(indexParameterEpoch.begin(), indexParameterEpoch.end(), Gnss::ParameterIndex());
    indexParameterBias = Gnss::ParameterIndex();

    // ===========================================================

    UInt use = useable() && normalEquationInfo.estimateReceiver.at(idRecv());
    Parallel::reduceSum(use, 0, normalEquationInfo.comm);
    Parallel::broadCast(use, 0, normalEquationInfo.comm);
    if(!use)
      return;

    // ===========================================================

    // Epoch parameter
    // ---------------
    if((normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_RECEIVER_KINEMATICPOSITION) ||
       (normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_RECEIVER_CLOCK))
      for(UInt idEpoch : normalEquationInfo.idEpochs)
        if(useable(idEpoch))
        {
          // kinematic position + clock
          // --------------------------
          std::vector<ParameterName> parameterNames;
          parameterNames.push_back(ParameterName(name(), "position.x", "", timeCorrected(idEpoch)));
          parameterNames.push_back(ParameterName(name(), "position.y", "", timeCorrected(idEpoch)));
          parameterNames.push_back(ParameterName(name(), "position.z", "", timeCorrected(idEpoch)));
          parameterNames.push_back(ParameterName(name(), "clock",      "", timeCorrected(idEpoch)));
          indexParameterEpoch.at(idEpoch) = normalEquationInfo.parameterNamesEpochReceiver(idEpoch, idRecv(), parameterNames);
        }

    // ===========================================================

    // signal bias
    // -----------
    if((normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_RECEIVER_SIGNALBIAS) ||
       (normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_RECEIVER_TECBIAS))
    {
      std::vector<ParameterName> parameterNames;
      Receiver::gnss().signalBiasParameter(signalBias.type, Receiver::gnss().groupTypes(Receiver::gnss().typesRecvTrans, idRecv(), TRUE/*isReceiver*/),
                                          /*eliminateClock*/(normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_RECEIVER_CLOCK),
                                          isCodeBiasEstimated(normalEquationInfo),
                                          isPhaseBiasEstimated(normalEquationInfo),
                                          /*estimateTecBias*/(normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_RECEIVER_TECBIAS),
                                          Bias, parameterNames);
      if(parameterNames.size())
      {
        for(auto &parameterName : parameterNames)
          parameterName.object = name();
        indexParameterBias = normalEquationInfo.parameterNamesReceiver(idRecv(), parameterNames);
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationReceiverLeo::aprioriParameter(const Gnss::NormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const
{
  try
  {
    // Epoch parameters
    // ----------------
    for(UInt idEpoch=0; idEpoch<indexParameterEpoch.size(); idEpoch++)
      if(indexParameterEpoch.at(idEpoch))
      {
        // kinematic position
        // ------------------
        copy((pos0.at(idEpoch)+dpos.at(idEpoch)).vector(), x0.row(normalEquationInfo.index(indexParameterEpoch.at(idEpoch)), 3));

        // clock
        // -----
        x0(normalEquationInfo.index(indexParameterEpoch.at(idEpoch)) + 3, 0) = LIGHT_VELOCITY * clockError(idEpoch);
      } // for(idEpoch)

    // signal biases
    // -------------
    if(indexParameterBias)
    {
      Vector l = signalBias.bias;
      Matrix A = Bias;
      copy(leastSquares(A, l), x0.row(normalEquationInfo.index(indexParameterBias), Bias.columns()));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double GnssParametrizationReceiverLeo::updateParameter(const Gnss::NormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef /*Wz*/, Bool printStatistics)
{
  try
  {
    // init convergence tracking
    // -------------------------
    Double maxChangeClock = 0;
    Double maxChangePos   = 0;
    Double maxChangeBias  = 0;
    std::string infoMaxChangeClock, infoMaxChangePos, infoMaxChangeBias;

    // Epoch parameters
    // ----------------
    for(UInt idEpoch=0; idEpoch<indexParameterEpoch.size(); idEpoch++)
      if(indexParameterEpoch.at(idEpoch))
      {
        const UInt index = normalEquationInfo.index(indexParameterEpoch.at(idEpoch));

        // kinematic position
        // ------------------
        const Vector3d dp(x(index+0,0), x(index+1,0), x(index+2,0));
        dpos.at(idEpoch) += dp;
        if(dp.r() > maxChangePos)
        {
          maxChangePos = dp.r();
          infoMaxChangePos = "  "+name()+": "+timeCorrected(idEpoch).dateTimeStr()+" posChange      = "+(1e3*dp.r())%"%6.1f"s+" mm";
        }
        if(dp.r() > 10000e3)
        {
          logWarning<<"  "<<stationInfo.markerName<<": "<<timeCorrected(idEpoch).dateTimeStr()<<" posChange   = "<<dp.r()<<" m, Epoch removed!"<<Log::endl;
          disable(idEpoch);
        }

        // clock
        // -----
        const Double clk = x(index+3, 0);
        updateClockError(idEpoch, clk/LIGHT_VELOCITY);
        if(std::fabs(clk) > maxChangeClock)
        {
           maxChangeClock = std::fabs(clk);
           infoMaxChangeClock = "  "+name()+": "+timeCorrected(idEpoch).dateTimeStr()+" clockChange    = "+(1e3*clk)%"%6.1f mm"s;
        }
        if(std::fabs(clk) > 10000e3)
        {
          logWarning<<"  "<<stationInfo.markerName<<": "<<timeCorrected(idEpoch).dateTimeStr()<<" clockChange = "<<clk<<" m, Epoch removed!"<<Log::endl;
          disable(idEpoch);
        }
      }

    // ===========================================================

    // update signal biases
    // --------------------
    if(indexParameterBias)
    {
      const Vector dBias = Bias * x.row(normalEquationInfo.index(indexParameterBias), Bias.columns());
      for(UInt idType=0; idType<dBias.size(); idType++)
        signalBias.bias.at(idType) += dBias(idType);
      for(UInt idType=0; idType<dBias.size(); idType++)
        if(std::fabs(dBias(idType)) > maxChangeBias)
        {
          maxChangeBias = std::fabs(dBias(idType));
          infoMaxChangeBias = "  "+name()+": "+signalBias.type.at(idType).str()+"              biasChange     = "+(1e3*dBias(idType))%"%6.1f mm"s;
        }
    }

    // ===========================================================

    Gnss::checkMaxChange(maxChangePos,   infoMaxChangePos,   printStatistics, normalEquationInfo.comm);
    Gnss::checkMaxChange(maxChangeClock, infoMaxChangeClock, printStatistics, normalEquationInfo.comm);
    Gnss::checkMaxChange(maxChangeBias,  infoMaxChangeBias,  printStatistics, normalEquationInfo.comm);

    return std::max(maxChangeClock, std::max(maxChangePos, maxChangeBias));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationReceiverLeo::updateCovariance(const Gnss::NormalEquationInfo &normalEquationInfo, const MatrixDistributed &covariance)
{
  try
  {
    // 3x3 epoch covariance matrix
    // ---------------------------
    for(UInt idEpoch=0; idEpoch<indexParameterEpoch.size(); idEpoch++)
      if(indexParameterEpoch.at(idEpoch))
      {
        const UInt idBlock = normalEquationInfo.block(indexParameterEpoch.at(idEpoch));
        const UInt index   = normalEquationInfo.index(indexParameterEpoch.at(idEpoch)) - covariance.blockIndex(idBlock);
        cov.at(idEpoch).xx() = covariance.N(idBlock, idBlock)(index+0, index+0);
        cov.at(idEpoch).yy() = covariance.N(idBlock, idBlock)(index+1, index+1);
        cov.at(idEpoch).zz() = covariance.N(idBlock, idBlock)(index+2, index+2);
        cov.at(idEpoch).xy() = covariance.N(idBlock, idBlock)(index+0, index+1);
        cov.at(idEpoch).xz() = covariance.N(idBlock, idBlock)(index+0, index+2);
        cov.at(idEpoch).yz() = covariance.N(idBlock, idBlock)(index+1, index+2);
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationReceiverLeo::writeResults(const Gnss::NormalEquationInfo &normalEquationInfo, const std::string &suffix)
{
  try
  {
    if(!useable())
      return;

    VariableList fileNameVariableList;
    addTimeVariables(fileNameVariableList);
    evaluateTimeVariables(0, timeIntervalStart, timeIntervalEnd, fileNameVariableList);

    // write kinematic orbits
    // ----------------------
    if(!fileNameOutOrbit.empty())
    {
      logStatus<<"write orbit data <"<<fileNameOutOrbit(fileNameVariableList).appendBaseName(suffix)<<">"<<Log::endl;
      OrbitArc arc;
      for(UInt idEpoch : normalEquationInfo.idEpochs)
        if(useable(idEpoch))
        {
          OrbitEpoch epoch;
          epoch.time      = timeCorrected(idEpoch);
          epoch.position  = position(idEpoch) - displacement.at(idEpoch);
          epoch.velocity  = velocity(idEpoch);
          epoch.acceleration.x() = clockError(idEpoch)*LIGHT_VELOCITY;
          arc.push_back(epoch);
        }
      InstrumentFile::write(fileNameOutOrbit(fileNameVariableList).appendBaseName(suffix), arc);
    }

    // Estimated clocks
    // ----------------
    if(!fileNameOutClock.empty())
    {
      logStatus<<"write clock error <"<<fileNameOutClock(fileNameVariableList).appendBaseName(suffix)<<">"<<Log::endl;
      MiscValueArc arc;
      for(UInt idEpoch : normalEquationInfo.idEpochs)
        if(useable(idEpoch))
        {
          const Double clk = clockError(idEpoch);
          if(clk)
          {
            MiscValueEpoch epoch;
            epoch.time  = timeCorrected(idEpoch);
            epoch.value = clk;
            arc.push_back(epoch);
          }
        }

      if(arc.size())
        InstrumentFile::write(fileNameOutClock(fileNameVariableList).appendBaseName(suffix), arc);
    }

    // write epoch covariance
    // ----------------------
    if(!fileNameOutCovariance.empty())
    {
      logStatus<<"write epoch covariance data <"<<fileNameOutCovariance(fileNameVariableList).appendBaseName(suffix)<<">"<<Log::endl;
      Covariance3dArc arc;
      for(UInt idEpoch : normalEquationInfo.idEpochs)
        if(useable(idEpoch))
        {
          Covariance3dEpoch epoch;
          epoch.time       = timeCorrected(idEpoch);
          epoch.covariance = cov.at(idEpoch);
          arc.push_back(epoch);
        }
      InstrumentFile::write(fileNameOutCovariance(fileNameVariableList).appendBaseName(suffix), arc);
    }

    // write residuals
    // ---------------
    if(!fileNameResiduals.empty())
    {
      logStatus<<"write GNSS receiver residuals <"<<fileNameResiduals(fileNameVariableList).appendBaseName(suffix)<<">"<<Log::endl;
      InstrumentFile::write(fileNameResiduals(fileNameVariableList).appendBaseName(suffix), residuals(normalEquationInfo.idEpochs));
    }

    // Signal bias
    // -----------
    if(!fileNameOutSignal.empty())
    {
      logStatus<<"write signal bias to file <"<<fileNameOutSignal(fileNameVariableList).appendBaseName(suffix)<<">"<<Log::endl;
      GnssSignalBias signalBias = this->signalBias;
      for(UInt idType=0; idType<signalBias.type.size(); idType++)
        if(signalBias.type.at(idType) == GnssType::PHASE)
          signalBias.bias.at(idType) = std::remainder(signalBias.bias.at(idType), LIGHT_VELOCITY/signalBias.type.at(idType).frequency());
      writeFileGnssSignalBias(fileNameOutSignal(fileNameVariableList).appendBaseName(suffix), signalBias);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

Bool GnssParametrizationReceiverLeo::isEpochEstimable(Gnss::AnalysisType analysisType, UInt idEpoch) const
{
  try
  {
    if(!useable(idEpoch))
      return FALSE;

    std::vector<GnssType> types;
    UInt count = 0;
    for(UInt idTrans=0; idTrans<Receiver::gnss().transmitter.size(); idTrans++)
      if(Receiver::gnss().transmitter.at(idTrans)->useable(idEpoch) && observation(idTrans, idEpoch) && observation(idTrans, idEpoch)->observationList(analysisType, types))
        count++;
    return (count>4); // more than 4 observation for overdetermination
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool GnssParametrizationReceiverLeo::isDesignMatrixReceiver(const Gnss::NormalEquationInfo &/*normalEquationInfo*/, UInt /*idTrans*/, UInt idEpoch) const
{
  try
  {
    if(!useable(idEpoch))
      return FALSE;
    return (indexParameterEpoch.at(idEpoch) || indexParameterBias);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationReceiverLeo::designMatrixReceiver(const Gnss::NormalEquationInfo &/*normalEquationInfo*/, const Gnss::ObservationEquation &eqn, Gnss::DesignMatrix &A) const
{
  try
  {
    // kinematic position + clock
    // --------------------------
    if(indexParameterEpoch.at(eqn.idEpoch))
      copy(eqn.A.column(Gnss::ObservationEquation::idxPosRecv,4), A.column(indexParameterEpoch.at(eqn.idEpoch)));

    // signal biases
    // -------------
    if(indexParameterBias)
    {
      MatrixSlice Design(A.column(indexParameterBias));
      for(UInt idType=0; idType<eqn.types.size(); idType++)
      {
        const UInt idx = GnssType::index(signalBias.type, eqn.types.at(idType));
        if(idx != NULLINDEX)
          matMult(1., eqn.A.column(Gnss::ObservationEquation::idxUnit + idType), Bias.row(idx), Design);
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Transform3d GnssParametrizationReceiverLeo::local2antennaFrame(UInt idEpoch) const
{
  try
  {
    return stationInfo.antenna.at(stationInfo.findAntenna(timeCorrected(idEpoch))).local2antennaFrame;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::vector<GnssType> GnssParametrizationReceiverLeo::definedTypes(UInt idEpoch) const
{
  try
  {
    UInt idRecv = stationInfo.findReceiver(timeCorrected(idEpoch));
    if(idRecv != NULLINDEX && stationInfo.receiver.at(idRecv).receiverDef)
      return stationInfo.receiver.at(idRecv).receiverDef->types;

    return std::vector<GnssType>();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector GnssParametrizationReceiverLeo::antennaVariations(UInt idEpoch, Angle azimut, Angle elevation, const std::vector<GnssType> &types) const
{
  try
  {
    Vector corr(types.size());
    corr += stationInfo.antennaVariations(timeCorrected(idEpoch), azimut, elevation, types, noPatternFoundAction);
    corr += signalBias.compute(types); // Code/Phase biases
    return corr;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector GnssParametrizationReceiverLeo::accuracy(UInt idEpoch, Angle azimut, Angle elevation, const std::vector<GnssType> &type) const
{
  try
  {
    return stationInfo.accuracy(timeCorrected(idEpoch), azimut, elevation, type, noPatternFoundAction);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool GnssParametrizationReceiverLeo::isCodeBiasEstimated(const Gnss::NormalEquationInfo &normalEquationInfo) const
{
  try
  {
    return estimateCodeBias && (normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_RECEIVER_SIGNALBIAS);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool GnssParametrizationReceiverLeo::isPhaseBiasEstimated(const Gnss::NormalEquationInfo &normalEquationInfo) const
{
  try
  {
    return estimatePhaseBias && (normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_RECEIVER_SIGNALBIAS);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
