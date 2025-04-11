/***********************************************/
/**
* @file gnssReceiverGeneratorLowEarthOrbiter.cpp
*
* @brief GNSS for Low Earth Orbiter (LEO).
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2021-02-25
*
*/
/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "inputOutput/logging.h"
#include "files/fileGnssSignalBias.h"
#include "files/filePlatform.h"
#include "files/fileInstrument.h"
#include "files/fileMatrix.h"
#include "gnss/gnss.h"
#include "gnss/gnssReceiver.h"
#include "gnss/gnssTransmitter.h"
#include "gnss/gnssReceiverGenerator/gnssReceiverGenerator.h"
#include "gnss/gnssReceiverGenerator/gnssReceiverGeneratorLowEarthOrbiter.h"

/***********************************************/

GnssReceiverGeneratorLowEarthOrbiter::GnssReceiverGeneratorLowEarthOrbiter(Config &config)
{
  try
  {
    std::string choice;

    readConfig(config, "inputfileStationInfo",         fileNameStationInfo,     Config::MUSTSET,  "{groopsDataDir}/gnss/receiverLowEarthOrbiter/", "satellite metadata (antenna, receiver, ...)");
    readConfig(config, "inputfileAntennaDefinition",   fileNameAntennaDef,      Config::MUSTSET,  "{groopsDataDir}/gnss/receiverLowEarthOrbiter/", "antenna center offsets and variations");
    if(readConfigChoice(config, "noAntennaPatternFound", choice, Config::MUSTSET, "ignoreObservation", "what should happen if no antenna pattern is found for an observation"))
    {
      if(readConfigChoiceElement(config, "ignoreObservation",   choice, "ignore observation if no matching pattern is found"))
        noPatternFoundAction = GnssAntennaDefinition::NoPatternFoundAction::IGNORE_OBSERVATION;
      if(readConfigChoiceElement(config, "useNearestFrequency", choice, "use pattern of nearest frequency if no matching pattern is found"))
        noPatternFoundAction = GnssAntennaDefinition::NoPatternFoundAction::USE_NEAREST_FREQUENCY;
      if(readConfigChoiceElement(config, "throwException",      choice, "throw exception if no matching pattern is found"))
        noPatternFoundAction = GnssAntennaDefinition::NoPatternFoundAction::THROW_EXCEPTION;
      endChoice(config);
    }
    readConfig(config, "inputfileReceiverDefinition",  fileNameReceiverDef,     Config::OPTIONAL, "{groopsDataDir}/gnss/receiverLowEarthOrbiter/", "observed signal types");
    readConfig(config, "inputfileAccuracyDefinition",  fileNameAccuracyDef,     Config::MUSTSET,  "{groopsDataDir}/gnss/receiverLowEarthOrbiter/", "elevation and azimut dependent accuracy");
    readConfig(config, "inputfileObservations",        fileNameObs,             Config::OPTIONAL,  "gnssReceiver_{loopTime:%D}.dat", "");
    readConfig(config, "inputfileOrbit",               fileNameOrbit,           Config::MUSTSET,  "",     "approximate positions");
    readConfig(config, "inputfileStarCamera",          fileNameStarCamera,      Config::MUSTSET,  "",     "satellite attitude");
    readConfig(config, "sigmaFactorPhase",             exprSigmaPhase,          Config::OPTIONAL, "",     "PHASE: factor = f(FREQ, ELE, SNR, ROTI, dTEc, IONOINDEX)");
    readConfig(config, "sigmaFactorCode",              exprSigmaCode,           Config::OPTIONAL, "",     "CODE: factor = f(FREQ, ELE, SNR, ROTI, dTEc, IONOINDEX)");
    readConfig(config, "supportsIntegerAmbiguities",   integerAmbiguities,      Config::DEFAULT,  "1",    "receiver tracks full cycle integer ambiguities");
    readConfig(config, "wavelengthFactor",             wavelengthFactor,        Config::DEFAULT,  "1.",   "factor to account for half-wavelength observations (collected by codeless squaring techniques)");
    readConfig(config, "useType",                      useType,                 Config::OPTIONAL, "",     "only use observations that match any of these patterns");
    readConfig(config, "ignoreType",                   ignoreType,              Config::OPTIONAL, "",     "ignore observations that match any of these patterns");
    readConfig(config, "elevationCutOff",              elevationCutOff,         Config::DEFAULT,  "0",    "[degree] ignore observations below cutoff");
    readConfig(config, "minObsCountPerTrack",          minObsCountPerTrack,     Config::DEFAULT,  "20",   "tracks with less number of epochs with observations are dropped");
    if(readConfigSequence(config, "preprocessing", Config::MUSTSET, "", "settings for preprocessing of observations/stations"))
    {
      readConfig(config, "printStatistics",              printInfo,               Config::DEFAULT,  "1",    "print preprocesssing statistics for all receivers");
      readConfig(config, "huber",                        huber,                   Config::DEFAULT,  "2.5",  "residuals > huber*sigma0 are downweighted");
      readConfig(config, "huberPower",                   huberPower,              Config::DEFAULT,  "1.5",  "residuals > huber: sigma=(e/huber)^huberPower*sigma0");
      readConfig(config, "codeMaxPositionDiff",          codeMaxPosDiff,          Config::DEFAULT,  "100",  "[m] max. allowed position error by PPP code only clock error estimation");
      readConfig(config, "denoisingLambda",              denoisingLambda,         Config::DEFAULT,  "5",    "regularization parameter for total variation denoising used in cylce slip detection");
      readConfig(config, "tecWindowSize",                tecWindowSize,           Config::DEFAULT,  "15",   "(0 = disabled) window size for TEC smoothness evaluation used in cycle slip detection");
      readConfig(config, "tecSigmaFactor",               tecSigmaFactor,          Config::DEFAULT,  "3.5",  "factor applied to moving standard deviation used as threshold in TEC smoothness evaluation during cycle slip detection");
      readConfig(config, "outputfileTrackBefore",        fileNameTrackBefore,     Config::OPTIONAL, "",     "variables {station}, {prn}, {timeStart}, {timeEnd}, {types}, TEC and MW-like combinations in cycles for each track before cycle slip detection");
      readConfig(config, "outputfileTrackAfter",         fileNameTrackAfter,      Config::OPTIONAL, "",     "variables {station}, {prn}, {timeStart}, {timeEnd}, {types}, TEC and MW-like combinations in cycles for each track after cycle slip detection");
      endSequence(config);
    } // readConfigSequence(preprocessing)
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssReceiverGeneratorLowEarthOrbiter::init(std::vector<GnssType> simulationTypes, const std::vector<Time> &times, const Time &timeMargin,
                                                const std::vector<GnssTransmitterPtr> &transmitters, EarthRotationPtr earthRotation,
                                                Parallel::CommunicatorPtr comm, std::vector<GnssReceiverPtr> &receivers)
{
  try
  {
    logStatus<<"init Low Earth Orbiter"<<Log::endl;
    Bool isSimulation = simulationTypes.size();
    if(isSimulation && !fileNameObs.empty())
      logWarningOnce<<"ignoring observation file <"<<fileNameObs<<">"<<Log::endl;

    std::vector<GnssAntennaDefinitionPtr> antennaDefList;
    if(!fileNameAntennaDef.empty())
      readFileGnssAntennaDefinition(fileNameAntennaDef, antennaDefList);

    std::vector<GnssReceiverDefinitionPtr> receiverDefList;
    if(!fileNameReceiverDef.empty())
      readFileGnssReceiverDefinition(fileNameReceiverDef, receiverDefList);

    std::vector<GnssAntennaDefinitionPtr> accuracyDefList;
    if(!fileNameAccuracyDef.empty())
      readFileGnssAntennaDefinition(fileNameAccuracyDef, accuracyDefList);

    Platform platform;
    readFilePlatform(fileNameStationInfo, platform);
    platform.name = platform.markerName;
    platform.fillGnssAntennaDefinition(antennaDefList);
    platform.fillGnssReceiverDefinition(receiverDefList);
    platform.fillGnssAccuracyDefinition(accuracyDefList);

    // test completeness of antennas
    for(const auto &instrument : platform.equipments)
    {
      auto antenna = std::dynamic_pointer_cast<PlatformGnssAntenna>(instrument);
      if(antenna && antenna->timeEnd > times.front() && antenna->timeStart <= times.back() && (!antenna->antennaDef || !antenna->accuracyDef))
        logWarningOnce<<platform.markerName<<"."<<platform.markerNumber<<": No "<<(!antenna->antennaDef ? "antenna" : "accuracy")<<" definition found for "<<antenna->str()<<Log::endl;
    }

    recv = std::make_shared<GnssReceiver>(Parallel::isMaster(comm), FALSE/*isEarthFixed*/, platform,
                                          noPatternFoundAction, Vector(times.size(), TRUE)/*useableEpochs*/,
                                          integerAmbiguities, wavelengthFactor);
    receivers.push_back(recv);

    if(Parallel::isMaster(comm))
    {
      try
      {
        recv->times = times;
        recv->clk.resize(times.size(), 0);
        recv->pos.resize(times.size());
        recv->vel.resize(times.size());
        recv->offset.resize(times.size());
        recv->global2local.resize(times.size());
        recv->global2antenna.resize(times.size());

        OrbitArc      orbit      = InstrumentFile::read(fileNameOrbit);
        StarCameraArc starCamera = InstrumentFile::read(fileNameStarCamera);
        Arc::checkSynchronized({orbit, starCamera});

        UInt i=0;
        for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
        {
          while((i < orbit.size()) && (orbit.at(i).time < times.at(idEpoch)-timeMargin))
            i++;
          if((i >= orbit.size()) || (orbit.at(i).time > times.at(idEpoch)+timeMargin))
          {
            recv->disable(idEpoch, "due to missing orbit/attitude data");
            continue;
          }

          auto antenna = platform.findEquipment<PlatformGnssAntenna>(times.at(idEpoch));
          if(!antenna || !antenna->antennaDef || !antenna->accuracyDef)
          {
            recv->disable(idEpoch, "missing antenna/accuracy patterns");
            continue;
          }

          recv->pos.at(idEpoch)            = orbit.at(i).position;
          recv->vel.at(idEpoch)            = orbit.at(i).velocity;
          recv->global2local.at(idEpoch)   = inverse(localNorthEastUp(recv->pos.at(idEpoch), Ellipsoid()));
          recv->global2antenna.at(idEpoch) = antenna->local2antennaFrame * inverse(starCamera.at(i).rotary);
          recv->offset.at(idEpoch)         = recv->global2local.at(idEpoch).transform(starCamera.at(i).rotary.rotate(antenna->position - platform.referencePoint(times.at(idEpoch))));
        }
        recv->preprocessingInfo("init()");

        logStatus<<"read observations"<<Log::endl;
        auto rotationCrf2Trf = std::bind(&EarthRotation::rotaryMatrix, earthRotation, std::placeholders::_1);
        if(isSimulation)
        {
          recv->simulateZeroObservations(simulationTypes, transmitters, rotationCrf2Trf, elevationCutOff,
                                         useType, ignoreType, GnssObservation::RANGE | GnssObservation::PHASE);
        }
        else
        {
          recv->readObservations(fileNameObs, transmitters,  rotationCrf2Trf, timeMargin, elevationCutOff,
                                 useType, ignoreType, GnssObservation::RANGE | GnssObservation::PHASE);
          // =================================
          if(recv->name() == "SENTINEL6A") //HACK: replace L2LG -> L2WG
          {
            logWarning<<recv->name()<<": replace L2LG -> L2WG"<<Log::endl;
            for(UInt idEpoch=0; idEpoch<recv->idEpochSize(); idEpoch++)
              for(UInt idTrans=0; idTrans<recv->idTransmitterSize(idEpoch); idTrans++)
              {
                GnssObservation *obs = recv->observation(idTrans, idEpoch);
                if(obs)
                  for(UInt idType=0; idType<obs->size(); idType++)
                    if(obs->at(idType).type == GnssType::L2_G+GnssType::L)
                      obs->at(idType).type = GnssType::L2_G + GnssType::W + (obs->at(idType).type & GnssType::PRN);
              }
          }
          // =================================
        }
      }
      catch(std::exception &e)
      {
        recv->disable(e.what());
      }
    } // if(isMaster())
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssReceiverGeneratorLowEarthOrbiter::preprocessing(Gnss *gnss, Parallel::CommunicatorPtr comm)
{
  try
  {
    if(fileNameObs.empty())
      return;

    logStatus<<"init observations"<<Log::endl;
    if(recv->isMyRank())
    {
      try
      {
        recv->createTracks(gnss->transmitters, minObsCountPerTrack, {GnssType::L5_G});
        std::vector<Vector3d> posApriori = recv->pos;
        recv->pos = recv->estimateInitialClockErrorFromCodeObservations(gnss->transmitters, gnss->funcRotationCrf2Trf, gnss->funcReduceModels, huber, huberPower, TRUE/*estimateKinematicPosition*/);
        // observation equations based on positions from code observations
        GnssReceiver::ObservationEquationList eqn(*recv, gnss->transmitters, gnss->funcRotationCrf2Trf, gnss->funcReduceModels, GnssObservation::RANGE | GnssObservation::PHASE);
        recv->pos = std::move(posApriori); // restore apriori positions

        recv->disableEpochsWithGrossCodeObservationOutliers(eqn, codeMaxPosDiff, 0.5);
        recv->writeTracks(fileNameTrackBefore, eqn, {GnssType::L5_G});
        recv->cycleSlipsDetection(eqn, minObsCountPerTrack, denoisingLambda, tecWindowSize, tecSigmaFactor, {GnssType::L5_G});
        recv->trackOutlierDetection(eqn, {GnssType::L5_G}, huber, huberPower);
        recv->cycleSlipsRepairAtSameFrequency(eqn);
        recv->writeTracks(fileNameTrackAfter, eqn, {GnssType::L5_G});

        // apply factors for accuracies from expressions
        if(exprSigmaPhase || exprSigmaCode)
        {
          for(UInt idTrans=0; idTrans<gnss->transmitters.size(); idTrans++)
            for(UInt idEpoch=0; idEpoch<recv->idEpochSize(); idEpoch++)
            {
              GnssObservation *obs = recv->observation(idTrans, idEpoch);
              if(obs)
              {
                VariableList varList;
                varList.setVariable("ROTI", 0);
                const UInt idx = obs->index(GnssType::ROTI);
                if(idx != NULLINDEX)
                  varList.setVariable("ROTI", obs->at(idx).observation);

                for(UInt idType=0; idType<obs->size(); idType++)
                  if((obs->at(idType).type == GnssType::RANGE) || (obs->at(idType).type == GnssType::PHASE))
                  {
                    varList.setVariable("FREQ",  obs->at(idType).type.frequency() );
                    varList.setVariable("SNR", 0);
                    const UInt idx = obs->index(GnssType::SNR + (obs->at(idType).type & GnssType::FREQUENCY));
                    if(idx != NULLINDEX)
                      varList.setVariable("SNR", obs->at(idx).observation);

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
        recv->disable(e.what());
      }
    } // if(recv->isMyRank())

    printPreprocessingInfos("preprocessing statistics after each step", {recv}, !printInfo/*disabledOnly*/, comm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssReceiverGeneratorLowEarthOrbiter::simulation(NoiseGeneratorPtr noiseClock, NoiseGeneratorPtr noiseObs,
                                                      Gnss *gnss, Parallel::CommunicatorPtr comm)
{
  try
  {
    logStatus<<"simulate observations"<<Log::endl;
    if(recv->isMyRank())
    {
      try
      {
        recv->simulateObservations(noiseClock, noiseObs, gnss->transmitters,
                                   gnss->funcRotationCrf2Trf, gnss->funcReduceModels,
                                   minObsCountPerTrack, Angle(0), GnssObservation::RANGE | GnssObservation::PHASE);
      }
      catch(std::exception &e)
      {
        recv->disable(e.what());
      }
    } // if(isMaster())

    printPreprocessingInfos("simulation statistics after each step", {recv}, !printInfo/*disabledOnly*/, comm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
