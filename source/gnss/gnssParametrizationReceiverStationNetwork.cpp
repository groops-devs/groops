/***********************************************/
/**
* @file gnssParametrizationReceiverStationNetwork.cpp
*
* @brief GNSS ground station network.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2010-08-19
*
*/
/***********************************************/

#include "base/import.h"
#include "base/planets.h"
#include "parser/dataVariables.h"
#include "config/config.h"
#include "inputOutput/logging.h"
#include "files/fileMatrix.h"
#include "files/fileParameterName.h"
#include "files/fileStringTable.h"
#include "files/fileInstrument.h"
#include "files/fileGnssSignalBias.h"
#include "files/fileGnssStationInfo.h"
#include "files/fileGriddedData.h"
#include "classes/troposphere/troposphere.h"
#include "classes/gravityfield/gravityfield.h"
#include "classes/tides/tides.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "misc/varianceComponentEstimation.h"
#include "gnss/gnss.h"
#include "gnss/gnssObservation.h"
#include "gnss/gnssDesignMatrix.h"
#include "gnss/gnssTransmitter.h"
#include "gnss/gnssParametrizationEarthRotation.h"
#include "gnss/gnssParametrizationReceiverStationNetwork.h"

/***********************************************/

GnssParametrizationReceiverStationNetwork::GnssParametrizationReceiverStationNetwork(Config &config)
{
  try
  {
    maxStationCount = MAX_UINT;
    noPatternFoundAction = GnssAntennaDefinition::NoPatternFoundAction::IGNORE_OBSERVATION;
    FileName fileNameStationList, fileNameStationInfo;
    FileName fileNameAntennaDef, fileNameReceiverDef, fileNameAccuracyDef;
    FileName deformationName, potentialName;
    Double   huber = 2.5, huberPower = 1.5;

    readConfig(config, "outputfileGriddedPosition",     fileNameGrid,                Config::OPTIONAL, "output/gridPosition_{loopTime:%D}.dat",              "delta north east up for all stations");
    readConfig(config, "outputfilePosition",            fileNamePosition,            Config::OPTIONAL, "output/stationPosition_{loopTime:%D}.{station}.dat", "full estimated coordinates (in TRF)");
    readConfig(config, "outputfilePositionSeries",      fileNamePositionSeries,      Config::OPTIONAL, "output/positionSeries_{loopTime:%D}.{station}.dat",  "time series of position deltas etc.");
    readConfig(config, "outputfileResiduals",           fileNameResiduals,           Config::OPTIONAL, "output/residuals_{loopTime:%D}.{station}.dat", "");
    readConfig(config, "outputfileSignalBias",          fileNameOutSignal,           Config::OPTIONAL, "output/signalBias_{loopTime:%D}.{station}.txt", "estimated signal biases");
    readConfig(config, "outputfileClock",               fileNameOutClock,            Config::OPTIONAL, "output/clock_{loopTime:%D}.{station}.dat", "estimated clock errors");
    readConfig(config, "outputfileUsedStationList",     fileNameUsedStationList,     Config::OPTIONAL, "", "ascii file with names of used stations");
    readConfig(config, "outputfileTroposphere",         fileNameOutTropo,            Config::OPTIONAL, "output/troposphere_{loopTime:%D}.{station}.txt", "columns: MJD, ZHD, ZWD, dry north gradient, wet north gradient, dry east gradient, wet east gradient");
    readConfig(config, "inputfileStationList",          fileNameStationList,         Config::MUSTSET,  "", "ascii file with station names");
    readConfig(config, "inputfileStationInfo",          fileNameStationInfo,         Config::MUSTSET,  "{groopsDataDir}/gnss/receiverStation/stationInfo/igs/stationInfo.{station}.xml", "station metadata (antennas, receivers, ...)");
    readConfig(config, "inputfileAntennaDefinition",    fileNameAntennaDef,          Config::MUSTSET,  "{groopsDataDir}/gnss/receiverStation/antennaDefinition/igs/igs14/antennaDefinition_igs14.dat", "antenna center offsets and variations");
    readConfig(config, "inputfileReceiverDefinition",   fileNameReceiverDef,         Config::OPTIONAL, "{groopsDataDir}/gnss/receiverStation/receiverDefinition/receiverDefinition.dat", "observed signal types");
    readConfig(config, "inputfileAccuracyDefinition",   fileNameAccuracyDef,         Config::MUSTSET,  "{groopsDataDir}/gnss/receiverStation/accuracyDefinition/accuracyDefinition.xml", "elevation and azimuth dependent accuracy");
    readConfig(config, "inputfileStationPosition",      fileNameStationPosition,     Config::OPTIONAL, "{groopsDataDir}/gnss/receiverStation/position/igs/igb14/stationPosition.{station}.dat", "precise coordinates used for no-net constraints (in TRF)");
    readConfig(config, "inputfileSignalBias",           fileNameSignalBias,          Config::OPTIONAL, "", "a priori signal biases");
    readConfig(config, "inputfileObservations",         fileNameObs,                 Config::OPTIONAL, "gnssReceiver_{loopTime:%D}.{station}.dat", "");
    readConfig(config, "kinematicPositionEstimation",   isKinematicPosition,         Config::DEFAULT,  "0", "positions are estimated every epoch (in TRF)");
    readConfig(config, "positionEstimation",            parametrizationPosition,     Config::DEFAULT,  "",  "parametrization of estimated position (in TRF)");
    readConfig(config, "noNetTranslationSigma",         sigmaNoNetTranslation,       Config::DEFAULT,  "0", "(0 = unconstrained) sigma [m] for no-net translation constraint on station coordinates");
    readConfig(config, "noNetRotationSigma",            sigmaNoNetRotation,          Config::DEFAULT,  "0", "(0 = unconstrained) sigma [m] at Earth's surface for no-net rotation constraint on station coordinates");
    readConfig(config, "noNetStationList",              fileNameNoNetStationList,    Config::OPTIONAL, "",  "ascii file with names of stations used for no-net constraints (empty = all stations)");
    readConfig(config, "estimateCodeBias",              estimateCodeBias,            Config::DEFAULT,  "0", "");
    readConfig(config, "estimatePhaseBias",             estimatePhaseBias,           Config::DEFAULT,  "1", "");
    readConfig(config, "supportsIntegerAmbiguities",    integerAmbiguities,          Config::DEFAULT,  "1", "receiver tracks full cycle integer ambiguities");
    readConfig(config, "troposphere",                   troposphere,                 Config::OPTIONAL, "",  "a priori troposphere model");
    readConfig(config, "troposphereWetEstimation",      parametrizationTropoWet,     Config::DEFAULT,  "",  "[m] parametrization of zenith wet delays");
    readConfig(config, "troposphereGradientEstimation", parametrizationTropoGradient,Config::DEFAULT,  "",  "[degree] parametrization of north and east gradients");
    readConfig(config, "antennaCenter",                 antennaCenterVariations,     Config::OPTIONAL, "",  "estimate antenna center variations");
    readConfig(config, "deformation",                   gravityfield,                Config::DEFAULT,  "",  "loading deformation");
    readConfig(config, "tides",                         tides,                       Config::DEFAULT,  "",  "tidal deformation");
    readConfig(config, "ephemerides",                   ephemerides,                 Config::OPTIONAL, "jpl", "for tidal deformation");
    readConfig(config, "inputfileDeformationLoadLoveNumber", deformationName,        Config::MUSTSET,  "{groopsDataDir}/loading/deformationLoveNumbers_CM_Gegout97.txt", "");
    readConfig(config, "inputfilePotentialLoadLoveNumber",   potentialName,          Config::OPTIONAL, "{groopsDataDir}/loading/loadLoveNumbers_Gegout97.txt", "if full potential is given and not only loading potential");
    if(readConfigSequence(config, "preprocessing", Config::MUSTSET, "", "settings for preprocessing of observations/stations"))
    {
      readConfig(config, "outputfileTrackBefore", fileNameTrackBefore,   Config::OPTIONAL, "",   "TEC and MW-like combinations in cycles for track before preprocessing (variables: loopTime, station, prn, trackTimeStart, trackTimeEnd)");
      readConfig(config, "outputfileTrackAfter",  fileNameTrackAfter,    Config::OPTIONAL, "",   "TEC and MW-like combinations in cycles for track after preprocessing (variables: loopTime, station, prn, trackTimeStart, trackTimeEnd)");
      readConfig(config, "maxStationCount",       maxStationCount,       Config::OPTIONAL, "",   "maximum number of stations to be used");
      readConfig(config, "elevationCutOff",       elevationCutOff,       Config::DEFAULT,  "5",  "[degree] ignore observations below cutoff");
      readConfig(config, "elevationTrackMinimum", elevationTrackMinimum, Config::DEFAULT,  "15", "[degree] ignore tracks that never exceed minimum elevation");
      readConfig(config, "useType",               useType,               Config::OPTIONAL, "",   "only use observations that match any of these patterns");
      readConfig(config, "ignoreType",            ignoreType,            Config::OPTIONAL, "",   "ignore observations that match any of these patterns");
      std::string choice;
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
      readConfig(config, "huber",                   huber,                   Config::DEFAULT,  "2.5",  "residuals > huber*sigma0 are downweighted");
      readConfig(config, "huberPower",              huberPower,              Config::DEFAULT,  "1.5",  "residuals > huber: sigma=(e/huber)^huberPower*sigma0");
      readConfig(config, "codeMaxPositionDiff",     codeMaxPosDiff,          Config::DEFAULT,  "100",  "[m] max. allowed position error by PPP code only clock error estimation");
      readConfig(config, "minObsCountPerTrack",     minObsCountPerTrack,     Config::DEFAULT,  "30",   "tracks with less number of epochs with observations are dropped");
      readConfig(config, "minEstimableEpochsRatio", minEstimableEpochsRatio, Config::DEFAULT,  "0.75", "[0,1] drop stations with lower ratio of estimable epochs to total epochs");
      readConfig(config, "denoisingLambda",         denoisingLambda,         Config::DEFAULT,  "5",    "regularization parameter for total variation denoising used in cylce slip detection");
      readConfig(config, "tecWindowSize",           tecWindowSize,           Config::DEFAULT,  "15",   "(0 = disabled) window size for TEC smoothness evaluation used in cycle slip detection");
      readConfig(config, "tecSigmaFactor",          tecSigmaFactor,          Config::DEFAULT,  "3.5",  "factor applied to moving standard deviation used as threshold in TEC smoothness evaluation during cycle slip detection");
      endSequence(config);
    } // readConfigSequence(preprocessing)
    if(isCreateSchema(config)) return;

    // ===========================================================

    std::vector<GnssAntennaDefinitionPtr> antennaDefList;
    if(!fileNameAntennaDef.empty())
      readFileGnssAntennaDefinition(fileNameAntennaDef, antennaDefList);

    std::vector<GnssReceiverDefinitionPtr> receiverDefList;
    if(!fileNameReceiverDef.empty())
      readFileGnssReceiverDefinition(fileNameReceiverDef, receiverDefList);

    std::vector<GnssAntennaDefinitionPtr> accuracyDefList;
    if(!fileNameAccuracyDef.empty())
      readFileGnssAntennaDefinition(fileNameAccuracyDef, accuracyDefList);

    // ===========================================================

    // init receivers
    // --------------
    VariableList fileNameVariableList;
    addVariable("station", fileNameVariableList);

    std::vector<std::vector<std::string>> stationName;
    readFileStringTable(fileNameStationList(fileNameVariableList), stationName);
    for(UInt i=0; i<stationName.size(); i++)
      for(UInt k=0; k<stationName.at(i).size(); k++) // alternatives
      {
        auto recv = std::make_shared<Receiver>();
        receiver.push_back(recv);

        // station info
        fileNameVariableList["station"]->setValue(stationName.at(i).at(k));
        readFileGnssStationInfo(fileNameStationInfo(fileNameVariableList), recv->stationInfo);

        recv->base = this;
        recv->stationName       = stationName.at(i).at(k);
        recv->countAlternatives = stationName.at(i).size()-1-k;
        recv->stationInfo.fillAntennaPattern (antennaDefList);
        recv->stationInfo.fillReceiverDefinition(receiverDefList);
        recv->stationInfo.fillAntennaAccuracy(accuracyDefList);
        recv->approxPosition = recv->stationInfo.approxPosition;
        if(recv->approxPosition.r()<6300e3)
          throw(Exception(recv->stationName+": "+recv->stationInfo.markerName+"."+recv->stationInfo.markerNumber+": No approx. position given"));
        recv->huber      = huber;
        recv->huberPower = huberPower;

        // local ellipsoidal north oriented frame (north, east, up)
        recv->lnof2trf = localNorthEastUp(recv->approxPosition, Ellipsoid());
      }

    // ===========================================================

    // load love numbers
    // -----------------
    if(!deformationName.empty())
    {
      Matrix love;
      readFileMatrix(deformationName, love);
      hn = love.column(0);
      ln = love.column(1);

      // models contain the total mass (loading mass & deformation mass effect)
      if(!potentialName.empty())
      {
        Vector kn;
        readFileMatrix(potentialName, kn);
        for(UInt n=2; n<std::min(kn.rows(), hn.rows()); n++)
          hn(n) /= (1.+kn(n));
        for(UInt n=2; n<std::min(kn.rows(), ln.rows()); n++)
          ln(n) /= (1.+kn(n));
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::vector<Gnss::ReceiverPtr> GnssParametrizationReceiverStationNetwork::receivers()
{
  std::vector<Gnss::ReceiverPtr> r;
  r.insert(r.begin(), receiver.begin(), receiver.end());
  return r;
}

/***********************************************/

std::vector<Gnss::ParametrizationPtr> GnssParametrizationReceiverStationNetwork::parametrizations()
{
  std::vector<Gnss::ParametrizationPtr> p({shared_from_this()});
  p.insert(p.begin(), antennaCenterVariations.begin(), antennaCenterVariations.end());
  return p;
}

/***********************************************/

void GnssParametrizationReceiverStationNetwork::initIntervalReceiver(Gnss::AnalysisType analysisType, const std::vector<Time> &times, const Time &timeMargin, Parallel::CommunicatorPtr comm)
{
  try
  {
    this->times = times;
    for(auto &recv : receiver)
      recv->indexParameterEpoch.resize(times.size(), Gnss::ParameterIndex());

    VariableList fileNameVariableList;
    addVariable("station", fileNameVariableList);

    // init stations (distributed)
    // ---------------------------
    // assign rank to stations with observations
    Vector receiverRank(receiver.size(), NULLINDEX);
    if(Parallel::isMaster(comm) && !fileNameObs.empty())
    {
      UInt count = 0;
      for(UInt idRecv=0; idRecv<receiver.size(); idRecv++)
      {
        try
        {
          fileNameVariableList["station"]->setValue(receiver.at(idRecv)->stationName);
          InstrumentFile(fileNameObs(fileNameVariableList));
          count++;
          receiverRank(idRecv) = 0;
          if(count > (maxStationCount+maxStationCount/10))
            break;
        }
        catch(...)
        {
        }
      } // for(idRecv)
      if(!count)
        throw(Exception("observation files <"+fileNameObs.str()+"> not found for any receiver"));
    }
    Parallel::broadCast(receiverRank, 0, comm);
    UInt idProcess = 0;
    for(UInt idRecv=0; idRecv<receiver.size(); idRecv++)
      if(receiverRank(idRecv) != NULLINDEX)
        receiverRank(idRecv) = (idProcess++) % Parallel::size(comm);

    // init stations only at one node
    for(UInt idRecv=0; idRecv<receiver.size(); idRecv++)
      if(Parallel::myRank(comm) == receiverRank(idRecv)) // only at one node
      {
        receiver.at(idRecv)->initInterval(times);
        receiver.at(idRecv)->posDisplacement.resize(times.size());
        receiver.at(idRecv)->dpos.resize(times.size());
        receiver.at(idRecv)->zenitDelayWet.resize(times.size(), 0.0);
        receiver.at(idRecv)->gradientX.resize(times.size(), 0.0);
        receiver.at(idRecv)->gradientY.resize(times.size(), 0.0);

        // test completeness of antennas
        for(const auto &antenna : receiver.at(idRecv)->stationInfo.antenna)
          if(antenna.timeEnd >= times.front() && antenna.timeStart < times.back() && (!antenna.antennaDef || !antenna.accuracyDef))
            logWarning<<receiver.at(idRecv)->stationInfo.markerName<<"."<<receiver.at(idRecv)->stationInfo.markerNumber<<":No "<<(!antenna.antennaDef ? "antenna" : "accuracy")
                      <<" definition found for "<<antenna.str()<<Log::endl;

        // disable epochs with incomplete antennas and calculate antenna center in TRF
        for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
        {
          const UInt idAnt = receiver.at(idRecv)->stationInfo.findAntenna(times.at(idEpoch));
          if(idAnt != NULLINDEX && receiver.at(idRecv)->stationInfo.antenna.at(idAnt).antennaDef && receiver.at(idRecv)->stationInfo.antenna.at(idAnt).accuracyDef)
            receiver.at(idRecv)->posDisplacement.at(idEpoch) += receiver.at(idRecv)->lnof2trf.transform(receiver.at(idRecv)->stationInfo.antenna.at(idAnt).position - receiver.at(idRecv)->stationInfo.referencePoint(times.at(idEpoch)));
          else
            receiver.at(idRecv)->disable(idEpoch);
        }
      } // for(idRecv)

    // ===========================================================

    // (optional) approximate station position
    // ---------------------------------------
    noNetStations = Vector(receiver.size(), TRUE);
    if(!fileNameStationPosition.empty())
    {
      for(UInt idRecv=0; idRecv<receiver.size(); idRecv++)
      {
        try
        {
          fileNameVariableList["station"]->setValue(receiver.at(idRecv)->stationName);
          Vector3dArc arc = InstrumentFile::read(fileNameStationPosition(fileNameVariableList));
          auto iter = arc.size() == 1 ? arc.begin() : std::find_if(arc.begin(), arc.end(), [&](const Epoch &e){ return e.time.isInInterval(times.front(), times.back()); });
          if(iter != arc.end())
            receiver.at(idRecv)->approxPosition = iter->vector3d;
          else
            noNetStations(idRecv) = FALSE;
        }
        catch(std::exception &/*e*/)
        {
          noNetStations(idRecv) = FALSE;
        }
      }
    }

    // ===========================================================

    // position list (TRF)
    // -------------------
    std::vector<Vector3d> positions;
    UInt idTropo = 0;
    for(auto &recv : receiver)
      if(recv->useable()) // use all for simulation
      {
        recv->idTropo = idTropo++;
        positions.push_back(recv->approxPosition);
      }

    // troposphere
    // -----------
    if(troposphere)
    {
      logStatus<<"init troposphere"<<Log::endl;
      troposphere->init(positions);
    }

    // tides & loading
    // ---------------
    if(positions.size())
    {
      logStatus<<"compute tides & loading"<<Log::endl;

      Vector gravity(positions.size()); // normal gravity
      for(UInt i=0; i<gravity.size(); i++)
        gravity(i) = Planets::normalGravity(positions.at(i));

      std::vector< std::vector<Vector3d> > disp(positions.size());
      for(UInt i=0; i<positions.size(); i++)
        disp.at(i).resize(times.size());

      std::vector<Rotary3d> rotEarth(times.size());
      for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
        rotEarth.at(idEpoch) = gnss().earthRotation->rotaryMatrix(times.at(idEpoch));

      tides->deformation(times, positions, rotEarth, gnss().earthRotation->earthRotation(), ephemerides, gravity, hn, ln, disp);
      gravityfield->deformation(times, positions, gravity, hn, ln, disp);
      tides        = nullptr;
      gravityfield = nullptr;

      // add displacements
      // -----------------
      for(auto &recv : receiver)
        if(recv->useable())
          for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
            recv->posDisplacement.at(idEpoch) += disp.at(recv->idTropo).at(idEpoch);
    } // if(positions.size())

    // ===========================================================

    // code & phase biases
    // -------------------
    if(!fileNameSignalBias.empty())
    {
      if(Parallel::isMaster(comm))
        logWarning<<"fileNameSignalBias not implemented"<<Log::endl;
//       for(auto &recv : receiver)
//         if(recv->useable())
//         {
//           fileNameVariableList["station"]->setValue(recv->stationName);
//           readFileGnssSignalBias(fileNameSignalBias(fileNameVariableList), recv->signalBias);
//         }
    }

    // ===========================================================

    // read observations
    // -----------------
    if(!fileNameObs.empty())
    {
      logStatus<<"read observations"<<Log::endl;
      logTimerStart;
      for(UInt idRecv=0; idRecv<receiver.size(); idRecv++)
        if(receiver.at(idRecv)->useable())
        {
          logTimerLoop(idRecv, receiver.size());
          fileNameVariableList["station"]->setValue(receiver.at(idRecv)->stationName);
          InstrumentFile fileReceiver(fileNameObs(fileNameVariableList));
          receiver.at(idRecv)->readObservations(fileReceiver, useType, ignoreType, timeMargin);
        } // for(idRecv)
      Parallel::barrier(comm);
      logTimerLoopEnd(receiver.size());
    }

    // ===========================================================

    // init observations
    // -----------------
    if(!fileNameObs.empty())
    {
      logStatus<<"init observations"<<Log::endl;
      Vector usedStations(receiver.size());
      Vector disabledStations(receiver.size());

      logTimerStart;
      for(UInt idRecv=0; idRecv<receiver.size(); idRecv++)
        if(receiver.at(idRecv)->useable())
        {
          logTimerLoop(idRecv, receiver.size());

          try
          {
            receiver.at(idRecv)->initObservation(analysisType);
            receiver.at(idRecv)->deleteUndefinedObservations();
            if(analysisType & Gnss::ANALYSIS_PHASE)
              receiver.at(idRecv)->createTracks(minObsCountPerTrack, {GnssType::L5_G});
            receiver.at(idRecv)->deleteObservationsOfInestimableEpochs(analysisType);
            receiver.at(idRecv)->estimateInitialClockErrorFromCodeObservations(codeMaxPosDiff, isKinematicPosition);

            Receiver::ObservationEquationList eqn(gnss(), *receiver.at(idRecv), analysisType & (Gnss::ANALYSIS_CODE | Gnss::ANALYSIS_PHASE));
            receiver.at(idRecv)->disableEpochsWithGrossCodeObservationOutliers(eqn, codeMaxPosDiff);
            if(analysisType & Gnss::ANALYSIS_PHASE)
            {
              receiver.at(idRecv)->removeTracksWithInsufficientObservationCount(eqn, minObsCountPerTrack, analysisType);
              fileNameVariableList["station"]->setValue(receiver.at(idRecv)->stationName);
              if(!fileNameTrackBefore.empty())
                receiver.at(idRecv)->writeTracks(fileNameTrackBefore, fileNameVariableList, eqn);
              receiver.at(idRecv)->cycleSlipsDetection(eqn, denoisingLambda, tecWindowSize, tecSigmaFactor);
              receiver.at(idRecv)->removeLowElevationTracks(eqn, elevationTrackMinimum);
              receiver.at(idRecv)->removeTracksWithInsufficientObservationCount(eqn, minObsCountPerTrack, analysisType);
              receiver.at(idRecv)->trackOutlierDetection(eqn, {GnssType::L5_G});
              receiver.at(idRecv)->cycleSlipsRepairAtSameFrequency(eqn);
              if(!fileNameTrackAfter.empty())
                receiver.at(idRecv)->writeTracks(fileNameTrackAfter, fileNameVariableList, eqn);
            }
          }
          catch(std::exception &e)
          {
            logWarning<<receiver.at(idRecv)->name()<<" disabled: "<<e.what()<<Log::endl;
            disabledStations(idRecv) = TRUE;
            receiver.at(idRecv)->disable();
          }

          // test estimable epochs
          if(!receiver.at(idRecv)->isReceiverEstimable(minEstimableEpochsRatio))
          {
            disabledStations(idRecv) = TRUE;
            receiver.at(idRecv)->disable();
          }

          usedStations(idRecv) = receiver.at(idRecv)->useable();
          if(usedStations(idRecv))
            for(UInt i=1; i<=receiver.at(idRecv)->countAlternatives; i++)
              receiver.at(idRecv+i)->disable();
        } // for(idRecv)
      Parallel::barrier(comm);
      logTimerLoopEnd(receiver.size());

      // disable alternative stations and limit to maxStationCount
      // ---------------------------------------------------------
      Parallel::reduceSum(usedStations, 0, comm);
      Parallel::broadCast(usedStations, 0, comm);
      UInt stationCount = 0;
      for(UInt idRecv=0; idRecv<receiver.size(); idRecv++)
      {
        if(usedStations(idRecv))
        {
          stationCount++;
          for(UInt i=1; i<=receiver.at(idRecv)->countAlternatives; i++)
          {
            receiver.at(idRecv+i)->disable();
            usedStations(idRecv+i)     = FALSE;
            disabledStations(idRecv+i) = FALSE;
          }
        }
        if(stationCount > maxStationCount)
        {
          receiver.at(idRecv)->disable();
          usedStations(idRecv)     = FALSE;
          disabledStations(idRecv) = FALSE;
        }
      }

      // Info about used stations
      // ------------------------
      Parallel::reduceSum(disabledStations, 0, comm);
      if(Parallel::isMaster(comm))
      {
        logInfo<<"  "<<sum(usedStations)<<" of "<<receiver.size()<<" stations used"<<Log::endl;
        if(sum(disabledStations))
        {
          std::vector<std::string> names;
          for(UInt idRecv=0; idRecv<receiver.size(); idRecv++)
            if(disabledStations(idRecv))
              names.push_back(receiver.at(idRecv)->name());
          std::sort(names.begin(), names.end());
          std::stringstream ss;
          for(const auto &name : names)
            ss<<name<<" ";
          logInfo<<"  "<<names.size()<<" disabled stations: "<<ss.str()<<Log::endl;
        }
      }

      // ===========================================================

      // calculate contribution of stations to net translation/rotation
      // --------------------------------------------------------------
      // info
      if(Parallel::isMaster(comm))
      {
        std::vector<std::string> names;
        for(UInt idRecv=0; idRecv<receiver.size(); idRecv++)
          if(usedStations(idRecv) && !noNetStations(idRecv))
            names.push_back(receiver.at(idRecv)->name());
        if(names.size())
        {
          std::stringstream ss;
          for(const auto &name : names)
            ss<<name<<" ";
          logInfo<<"  "<<names.size()<<" stations without precise position: "<<ss.str()<<Log::endl;
          if(sigmaNoNetRotation || sigmaNoNetTranslation)
            logInfo<<"these stations are excluded from no-net rotation/translation constraints"<< Log::endl;
        }
      }

      // exclude unuseable stations
      for(UInt idRecv=0; idRecv<receiver.size(); idRecv++)
        if(!usedStations(idRecv))
          noNetStations(idRecv) = FALSE;

      if(!fileNameNoNetStationList.empty())
      {
        Vector noNetStationsAll = noNetStations;
        noNetStations = Vector(receiver.size());

        std::vector<std::vector<std::string>> noNetStationName;
        readFileStringTable(fileNameNoNetStationList(fileNameVariableList), noNetStationName);
        for(UInt i=0; i<noNetStationName.size(); i++)
          for(UInt k=0; k<noNetStationName.at(i).size(); k++) // alternatives
          {
            const UInt idRecv = std::distance(receiver.begin(), std::find_if(receiver.begin(), receiver.end(), [&](const Gnss::ReceiverPtr &recv) {return (recv->name() == noNetStationName.at(i).at(k));}));
            if((idRecv < receiver.size()) && noNetStationsAll(idRecv))
            {
              noNetStations(idRecv) = TRUE;
              break;
            }
          }
      }

      netRotation    = Vector(3);
      netTranslation = Vector(3);
      if(sigmaNoNetRotation || sigmaNoNetTranslation)
      {
        if(sum(noNetStations) > 0)
          logInfo<<"  "<<sum(noNetStations)<<" stations contribute to the computation of net translation/rotation"<<Log::endl;
        else
          throw(Exception("no stations contribute to the computation of net translation/rotation"));
      }
    } // if(fileNameObs)

    // ===========================================================

    // init temporal parametrization
    // -----------------------------
    if(parametrizationPosition)
    {
      parametrizationPosition->setInterval(times.front(), times.back(), TRUE);
      const UInt countParameterPosition = 3*parametrizationPosition->parameterCount();
      for(UInt idRecv=0; idRecv<receiver.size(); idRecv++)
        receiver.at(idRecv)->xPos = Vector(countParameterPosition);
    }

    if(parametrizationTropoWet)
    {
      parametrizationTropoWet->setInterval(times.front(), times.back(), TRUE);
      const UInt countParameterTropoWet = parametrizationTropoWet->parameterCount();
      for(UInt idRecv=0; idRecv<receiver.size(); idRecv++)
        receiver.at(idRecv)->xTropoWet = Vector(countParameterTropoWet);
    }

    if(parametrizationTropoGradient)
    {
      parametrizationTropoGradient->setInterval(times.front(), times.back(), TRUE);
      const UInt countParameterTropoGradient = 2*parametrizationTropoGradient->parameterCount();
      for(UInt idRecv=0; idRecv<receiver.size(); idRecv++)
        receiver.at(idRecv)->xTropoGradient = Vector(countParameterTropoGradient);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool GnssParametrizationReceiverStationNetwork::initIntervalDisabling(Gnss::AnalysisType analysisType, const std::vector<Time> &/*times*/, const Time &/*timeMargin*/, Parallel::CommunicatorPtr comm)
{
  try
  {
    UInt disabled = FALSE;
    Gnss::Receiver::ObservationEquationList eqn;
    for(auto recv : receiver)
      if(recv->useable())
      {
        disabled = recv->deleteObservationsOfInestimableEpochs(analysisType) || disabled;
        disabled = recv->removeTracksWithInsufficientObservationCount(eqn, minObsCountPerTrack, analysisType) || disabled;

        // test estimable epochs
        if(!recv->isReceiverEstimable(minEstimableEpochsRatio))
        {
          logWarning<<recv->name()<<" disabled: not enough estimable epochs"<<Log::endl;
          recv->disable();
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

void GnssParametrizationReceiverStationNetwork::initIntervalLate(Gnss::AnalysisType /*analysisType*/, const std::vector<Time> &/*times*/, const Time &/*timeMargin*/, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    for(auto &recv : receiver)
    {
      // get observation types
      std::vector<GnssType> types;
      for(auto &typesTrans : gnss().typesRecvTrans.at(recv->idRecv()))
        for(GnssType type : typesTrans)
          if((type == GnssType::PHASE) || (type == GnssType::RANGE))
            if(GnssType::index(types, type) == NULLINDEX)
              types.push_back(type);
      if(types.size() == 0)
        continue;
      std::sort(types.begin(), types.end());

      // code & phase biases
      // -------------------
      if(types.size()) // consider simulation case
      {
        recv->signalBias.bias = recv->signalBias.compute(types); // apriori signal bias
        recv->signalBias.type = types;
      }

      // antenna center variations
      // -------------------------
      for(auto &acv : antennaCenterVariations)
        acv->addAntennas(recv->idRecv(), TRUE/*deviceIsReceiver*/, recv->stationInfo, types);
    } // for(idRecv)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationReceiverStationNetwork::initParameter(Gnss::NormalEquationInfo &normalEquationInfo)
{
  try
  {
    // reset parameter indices
    // -----------------------
    for(UInt idRecv=0; idRecv<receiver.size(); idRecv++)
    {
      std::fill(receiver.at(idRecv)->indexParameterEpoch.begin(), receiver.at(idRecv)->indexParameterEpoch.end(), Gnss::ParameterIndex());
      receiver.at(idRecv)->indexParameterPosition      = Gnss::ParameterIndex();
      receiver.at(idRecv)->indexParameterTropoWet      = Gnss::ParameterIndex();
      receiver.at(idRecv)->indexParameterTropoGradient = Gnss::ParameterIndex();
      receiver.at(idRecv)->indexParameterBias          = Gnss::ParameterIndex();
    }

    // ===========================================================

    // distribute usable epochs
    Matrix useableEpoch(receiver.size(), times.size());
    for(UInt idRecv=0; idRecv<receiver.size(); idRecv++)
      if(normalEquationInfo.estimateReceiver.at(receiver.at(idRecv)->idRecv()) && receiver.at(idRecv)->useable())
        for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
          useableEpoch(idRecv, idEpoch) = receiver.at(idRecv)->isEpochEstimable(normalEquationInfo.analysisType, idEpoch);
    Parallel::reduceSum(useableEpoch, 0, normalEquationInfo.comm);
    Parallel::broadCast(useableEpoch, 0, normalEquationInfo.comm);

    for(UInt idRecv=0; idRecv<receiver.size(); idRecv++)
      if(sum(useableEpoch.row(idRecv)))
      {
        // Epoch parameter
        // ---------------
        for(UInt idEpoch : normalEquationInfo.idEpochs)
          if(useableEpoch(idRecv, idEpoch))
          {
            std::vector<ParameterName> parameterNames;

            // clock
            // -----
            if(normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_RECEIVER_CLOCK)
              parameterNames.push_back(ParameterName(receiver.at(idRecv)->name(), "clock", "", times.at(idEpoch)));

            // kinematic position
            // ------------------
            if(isKinematicPosition && (normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_RECEIVER_KINEMATICPOSITION))
            {
              parameterNames.push_back(ParameterName(receiver.at(idRecv)->name(), "position.x", "", times.at(idEpoch)));
              parameterNames.push_back(ParameterName(receiver.at(idRecv)->name(), "position.y", "", times.at(idEpoch)));
              parameterNames.push_back(ParameterName(receiver.at(idRecv)->name(), "position.z", "", times.at(idEpoch)));
            }

            receiver.at(idRecv)->indexParameterEpoch.at(idEpoch) = normalEquationInfo.parameterNamesEpochReceiver(idEpoch, idRecv, parameterNames);
          } // for(idEpoch)


        // position
        // --------
        if(parametrizationPosition->parameterCount() && !isKinematicPosition && (normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_RECEIVER_POSITION))
        {
          std::vector<ParameterName> parameterNames;
          std::vector<ParameterName> name({{receiver.at(idRecv)->name(), "position.x"}, {receiver.at(idRecv)->name(), "position.y"}, {receiver.at(idRecv)->name(), "position.z"}});
          parametrizationPosition->parameterName(name, parameterNames);
          receiver.at(idRecv)->indexParameterPosition = normalEquationInfo.parameterNamesReceiver(receiver.at(idRecv)->idRecv(), parameterNames);
        }

        // wet troposphere
        // ---------------
        if(parametrizationTropoWet->parameterCount() && (normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_RECEIVER_TROPOSPHERE_WET))
        {
          std::vector<ParameterName> parameterNames;
          parametrizationTropoWet->parameterName({ParameterName(receiver.at(idRecv)->name(), "troposphereWet")}, parameterNames);
          receiver.at(idRecv)->indexParameterTropoWet = normalEquationInfo.parameterNamesReceiver(receiver.at(idRecv)->idRecv(), parameterNames);
        }

        // troposphere gradient
        // --------------------
        if(parametrizationTropoGradient->parameterCount() && (normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_RECEIVER_TROPOSPHERE_GRADIENT))
        {
          std::vector<ParameterName> parameterNames;
          std::vector<ParameterName> name({{receiver.at(idRecv)->name(), "troposphereGradient.x"}, {receiver.at(idRecv)->name(), "troposphereGradient.y"}});
          parametrizationTropoGradient->parameterName(name, parameterNames);
          receiver.at(idRecv)->indexParameterTropoGradient = normalEquationInfo.parameterNamesReceiver(receiver.at(idRecv)->idRecv(), parameterNames);
        }

        // signal bias
        // -----------
        if((normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_RECEIVER_SIGNALBIAS) ||
           (normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_RECEIVER_TECBIAS))
        {
          std::vector<ParameterName> parameterNames;
          gnss().signalBiasParameter(receiver.at(idRecv)->signalBias.type, gnss().groupTypes(gnss().typesRecvTrans, receiver.at(idRecv)->idRecv(), TRUE/*isReceiver*/),
                                      /*eliminateClock*/(normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_RECEIVER_CLOCK),
                                      receiver.at(idRecv)->isCodeBiasEstimated(normalEquationInfo),
                                      receiver.at(idRecv)->isPhaseBiasEstimated(normalEquationInfo),
                                      /*estimateTecBias*/(normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_RECEIVER_TECBIAS),
                                      receiver.at(idRecv)->Bias, parameterNames);
          if(parameterNames.size())
          {
            for(auto &parameterName : parameterNames)
              parameterName.object = receiver.at(idRecv)->name();
            receiver.at(idRecv)->indexParameterBias = normalEquationInfo.parameterNamesReceiver(receiver.at(idRecv)->idRecv(), parameterNames);
          }
        }
      } // for(idRecv)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationReceiverStationNetwork::aprioriParameter(const Gnss::NormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const
{
  try
  {
    for(auto &recv : receiver)
      if(recv->useable())
      {
        // Epoch parameters
        // ----------------
        for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
        {
          // clock
          // -----
          if(recv->indexParameterEpoch.at(idEpoch))
            x0(normalEquationInfo.index(recv->indexParameterEpoch.at(idEpoch)), 0) = LIGHT_VELOCITY * recv->clockError(idEpoch);

          // kinematic position
          // ------------------
          if(isKinematicPosition && (recv->indexParameterEpoch.at(idEpoch)))
            copy((recv->approxPosition+recv->dpos.at(idEpoch)).vector(), x0.row(normalEquationInfo.index(recv->indexParameterEpoch.at(idEpoch))+1, 3));
        }// for(idEpoch)

        // positions
        // ---------
        if(recv->indexParameterPosition)
        {
          Vector l(3*times.size());
          Matrix A(3*times.size(), 3*parametrizationPosition->parameterCount());
          for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
          {
            copy((recv->approxPosition+recv->dpos.at(idEpoch)).vector(), l.row(3*idEpoch, 3));
            parametrizationPosition->designMatrix(times.at(idEpoch), identityMatrix(3), A.row(3*idEpoch, 3));
          }
          copy(leastSquares(A, l), x0.row(normalEquationInfo.index(recv->indexParameterPosition), 3*parametrizationPosition->parameterCount()));
        }

        // wet troposphere
        // ---------------
        if(recv->indexParameterTropoWet)
          copy(recv->xTropoWet, x0.row(normalEquationInfo.index(recv->indexParameterTropoWet), parametrizationTropoWet->parameterCount()));

        // troposphere gradient
        // --------------------
        if(recv->indexParameterTropoGradient)
          axpy(RAD2DEG, recv->xTropoGradient, x0.row(normalEquationInfo.index(recv->indexParameterTropoGradient), 2*parametrizationTropoGradient->parameterCount()));

        // signal biases
        // -------------
        if(recv->indexParameterBias)
        {
          Vector l = recv->signalBias.bias;
          Matrix A = recv->Bias;
          copy(leastSquares(A, l), x0.row(normalEquationInfo.index(recv->indexParameterBias), recv->Bias.columns()));
        }
      } // for(idRecv)

    // ===========================================================
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationReceiverStationNetwork::observationEquation(const Gnss::NormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const
{
  try
  {
    // zero mean code bias for clock for each satellite group (relative to first one)
    // ------------------------------------------------------------------------------
    if(Parallel::isMaster(normalEquationInfo.comm) && (normalEquationInfo.estimationType & Gnss::NormalEquationInfo::CONSTRAINT_RECEIVER_SIGNALBIAS_DIFFCLOCK) &&
       std::any_of(receiver.begin(), receiver.end(), [](const auto &recv){return recv->indexParameterBias;}))
    {
      // for each receiver and transmitter: used types (receiver types)
      std::vector<std::vector<std::vector<GnssType>>> typesRecvTrans;
      for(auto &recv : receiver)
        if(recv->indexParameterBias)
          typesRecvTrans.push_back(gnss().typesRecvTrans.at(recv->idRecv()));

      std::vector<std::vector<GnssType>> groupTypes = gnss().groupTypes(typesRecvTrans, NULLINDEX/*idRecv*/, TRUE/*isReceiver*/);
      if(groupTypes.size() > 1)
      {
        // simulate simplified normal equations (each receiver observes satellites groups)
        // parameters: all receiver biases, receiver/transmitter clocks, STEC
        Matrix N;
        std::vector<std::vector<GnssType>> recvTypes;
        std::vector<UInt>                  indexBias;
        {
          UInt countBias = 0;
          std::vector<Matrix> Bias;
          for(auto &recv : receiver)
            if(recv->indexParameterBias)
            {
              if(find(recvTypes.begin(), recvTypes.end(), recv->signalBias.type) != recvTypes.end())
                continue;
              Bias.push_back(recv->Bias);
              recvTypes.push_back(recv->signalBias.type);
              indexBias.push_back(countBias);
              countBias += recv->Bias.columns();
            }
          UInt countRecv = recvTypes.size();

          Matrix N11(countRecv + groupTypes.size(), Matrix::SYMMETRIC); // receiver and transmitter clocks
          Matrix N22(countBias, Matrix::SYMMETRIC);                     // receiver biases
          Matrix N12(countRecv + groupTypes.size(), countBias);

          UInt idxRecv = 0;
          for(UInt idRecv=0; idRecv<recvTypes.size(); idRecv++)
          {
            UInt   idxBias = indexBias.at(idRecv);
            Matrix B = Bias.at(idRecv);
            Matrix ClkRecv(B.rows(), 1, 1.); // receiver clock
            Matrix ClkTrans(B.rows(), groupTypes.size());
            Matrix Tec(B.rows(), groupTypes.size());
            UInt countTec = 0;
            for(UInt idGroup=0; idGroup<groupTypes.size(); idGroup++)
            {
              Bool used = FALSE;
              for(GnssType type : groupTypes.at(idGroup))
              {
                const UInt idx = GnssType::index(recvTypes.at(idRecv), type);
                if(idx == NULLINDEX)
                  continue;
                used = TRUE;
                Tec(idx, countTec) = type.ionosphericFactor();
                ClkTrans(idx, idGroup) = 1.; // transmitter clock
              }
              if(used)
                countTec++;
            }
            // eliminate Tec
            eliminationParameter(Tec.column(0, countTec), {ClkRecv, ClkTrans, B});
            // accumulate normals
            rankKUpdate(1, ClkRecv,  N11.slice(idxRecv, idxRecv, ClkRecv.columns(), ClkRecv.columns()));
            rankKUpdate(1, ClkTrans, N11.slice(countRecv, countRecv, ClkTrans.columns(), ClkTrans.columns()));
            rankKUpdate(1, B,     N22.slice(idxBias, idxBias, B.columns(), B.columns()));
            matMult(1., ClkRecv.trans(),  ClkTrans, N11.slice(idxRecv,   countRecv, ClkRecv.columns(),  ClkTrans.columns()));
            matMult(1., ClkRecv.trans(),  B,     N12.slice(idxRecv,   idxBias,   ClkRecv.columns(),  B.columns()));
            matMult(1., ClkTrans.trans(), B,     N12.slice(countRecv, idxBias,   ClkTrans.columns(), B.columns()));
            idxRecv++;
          }
          // zero mean of transmitter clocks
          rankKUpdate(1, Vector(groupTypes.size(), 1.).trans(), N11.slice(countRecv, countRecv, groupTypes.size(), groupTypes.size()));
          // eliminate clock parameters
          cholesky(N11);
          triangularSolve(1., N11.trans(), N12);
          rankKUpdate(-1, N12, N22);
          // count zero eigen values => null space
          Vector eigen = eigenValueDecomposition(N22, TRUE/*computeEigenVectors*/);
          UInt count = 0;
          while((count<eigen.rows()) && (eigen(count) < 1e-5))
            count++;
          // null space defines the constraint equations
          N = N22.column(0, count).trans();
        } // simulation

        logStatus<<"apply "<<N.rows()<<" zero mean equations of receiver code biases to align "<<groupTypes.size()<<" different satellite groups"<<Log::endl;
        if(N.rows())
        {
          const Double sigma = 0.000001; // 0.1 mm
          Gnss::DesignMatrix A(normalEquationInfo,  Vector(N.rows()));
          for(auto &recv : receiver)
            if(recv->indexParameterBias)
            {
              const UInt idx = indexBias.at(std::distance(recvTypes.begin(), find(recvTypes.begin(), recvTypes.end(), recv->signalBias.type)));
              axpy(1./sigma, N.column(idx, recv->Bias.columns()), A.column(recv->indexParameterBias));
            }
          A.accumulateNormals(normals, n, lPl, obsCount);
        }
      } // if(groupTypes.size() > 1)
    }

    // ===========================================================

    // no-net constrains
    // -----------------
    const Bool applyNoNetTranslation = sigmaNoNetTranslation && (normalEquationInfo.estimationType & Gnss::NormalEquationInfo::CONSTRAINT_RECEIVER_POSITION_NONETTRANSLATION);
    const Bool applyNoNetRotation    = sigmaNoNetRotation    && (normalEquationInfo.estimationType & Gnss::NormalEquationInfo::CONSTRAINT_RECEIVER_POSITION_NONETROTATION);
    if(applyNoNetTranslation || applyNoNetRotation)
    {
      std::vector<UInt> noNetReceiverIndexes;
      for(UInt idRecv=0; idRecv<receiver.size(); idRecv++)
        if(noNetStations(idRecv) && receiver.at(idRecv)->indexParameterPosition)
          noNetReceiverIndexes.push_back(idRecv);

      // collect l vector for weight estimation
      Vector l(3*noNetReceiverIndexes.size());
      for(UInt i=0; i<noNetReceiverIndexes.size(); i++)
        copy(receiver.at(noNetReceiverIndexes.at(i))->xPos, l.row(3*i, 3));
      Parallel::reduceSum(l, 0, normalEquationInfo.comm);

      if(Parallel::isMaster(normalEquationInfo.comm))
      {
        if(applyNoNetTranslation) logStatus<<"apply no-net translation to receiver positions"<<Log::endl;
        if(applyNoNetRotation)    logStatus<<"apply no-net rotation to receiver positions"<<Log::endl;

        UInt obsCount = 0;
        UInt idxNNT = obsCount; if(applyNoNetTranslation) obsCount += 3;
        UInt idxNNR = obsCount; if(applyNoNetRotation)    obsCount += 3;

        // compute station weights for no-net constraint
        Matrix A(l.size(), obsCount);
        for(UInt i=0; i<noNetReceiverIndexes.size(); i++)
        {
          if(applyNoNetTranslation)
            copy(identityMatrix(3), A.slice(3*i, idxNNT, 3, 3));
          if(applyNoNetRotation)
          {
            const Vector3d pos = receiver.at(noNetReceiverIndexes.at(i))->approxPosition/DEFAULT_R;
            A(3*i+0, idxNNR+1) =  pos.z(); A(3*i+0, idxNNR+2) = -pos.y();
            A(3*i+1, idxNNR+0) = -pos.z(); A(3*i+1, idxNNR+2) =  pos.x();
            A(3*i+2, idxNNR+0) =  pos.y(); A(3*i+2, idxNNR+1) = -pos.x();
          }
        } // for(recv)
        Vector sigmas;
        Vce::robustLeastSquares(A, l, 3, receiver.at(noNetReceiverIndexes.front())->huber, receiver.at(noNetReceiverIndexes.front())->huberPower, 5, sigmas);
        Double totalWeight = 0;
        for(UInt i=0; i<noNetReceiverIndexes.size(); i++)
          totalWeight += 1./std::pow(sigmas(i), 2);

        // weighted no-net constraints
        Vector l2(obsCount);
        if(applyNoNetTranslation) axpy(-1./sigmaNoNetTranslation, netTranslation, l2.row(idxNNT, 3));  // constrain towards zero (0-x0)
        if(applyNoNetRotation)    axpy(-1./sigmaNoNetRotation,    netRotation,    l2.row(idxNNR, 3));  // constrain towards zero (0-x0)

        Gnss::DesignMatrix A2(normalEquationInfo, l2);
        for(UInt i=0; i<noNetReceiverIndexes.size(); i++)
        {
          Matrix T(obsCount, 3);
          if(applyNoNetTranslation)
            axpy(1./sigmaNoNetTranslation, identityMatrix(3), T.row(idxNNT, 3));
          if(applyNoNetRotation)
          {
            const Vector3d pos = receiver.at(noNetReceiverIndexes.at(i))->approxPosition * (1./sigmaNoNetRotation/DEFAULT_R);
            T(idxNNR+0, 1) =  pos.z(); T(idxNNR+0, 2) = -pos.y();
            T(idxNNR+1, 0) = -pos.z(); T(idxNNR+1, 2) =  pos.x();
            T(idxNNR+2, 0) =  pos.y(); T(idxNNR+2, 1) = -pos.x();
          }

          // apply no net for every temporal base function (constant, trend, ...)
          MatrixSlice Design(A2.column(receiver.at(noNetReceiverIndexes.at(i))->indexParameterPosition));
          for(UInt j=0; j<parametrizationPosition->parameterCount(); j++)
            axpy(1./std::pow(sigmas(i), 2)/totalWeight, T, Design.column(3*j, 3));
        } // for(recv)
        A2.accumulateNormals(normals, n, lPl, obsCount);
      }
    } // if(no-net constrains)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double GnssParametrizationReceiverStationNetwork::updateParameter(const Gnss::NormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef /*Wz*/, Bool printStatistics)
{
  try
  {
    // init convergence tracking
    // -------------------------
    Double      maxChangeClock    = 0;
    Double      maxChangePos      = 0;
    Double      maxChangeWet      = 0;
    Double      maxChangeGradient = 0;
    Double      maxChangeBias     = 0;
    std::string infoMaxChangeClock, infoMaxChangePos, infoMaxChangeWet, infoMaxChangeGradient, infoMaxChangeBias;

    for(auto &recv : receiver)
      if(recv->useable())
      {
        // Epoch parameters
        // ----------------
        for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
        {
          // clock
          // -----
          if(recv->indexParameterEpoch.at(idEpoch))
          {
            const Double clk = x(normalEquationInfo.index(recv->indexParameterEpoch.at(idEpoch)), 0);
            recv->updateClockError(idEpoch, clk/LIGHT_VELOCITY);
            if(std::fabs(clk) > maxChangeClock)
            {
              maxChangeClock = std::fabs(clk);
              infoMaxChangeClock = "  "+recv->name()+": "+recv->timeCorrected(idEpoch).dateTimeStr()+" clockChange    = "+(1e3*clk)%"%6.1f mm"s;
            }
          }

          // kinematic position
          // ------------------
          if(isKinematicPosition && recv->indexParameterEpoch.at(idEpoch))
          {
            const Vector3d dp = Vector3d(x.row(normalEquationInfo.index(recv->indexParameterEpoch.at(idEpoch))+1, 3));
            recv->dpos.at(idEpoch) += dp;
            if(dp.r() > maxChangePos)
            {
              maxChangePos = dp.r();
              infoMaxChangePos = "  "+recv->name()+": "+recv->timeCorrected(idEpoch).dateTimeStr()+" posChange      = "+(1e3*maxChangePos)%"%6.1f mm"s;
            }
          }
        }// for(idEpoch)

        // update positions
        // ----------------
        if(recv->indexParameterPosition)
        {
          recv->xPos += x.row(normalEquationInfo.index(recv->indexParameterPosition), 3*parametrizationPosition->parameterCount());
          std::vector<UInt>   index;
          std::vector<Double> factor;
          for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
          {
            parametrizationPosition->factors(std::max(recv->timeCorrected(idEpoch), times.at(0)), index, factor);
            Vector p(3);
            for(UInt k=0; k<factor.size(); k++)
              axpy(factor.at(k), recv->xPos.row(3*index.at(k),3), p);
            Vector3d posOld = recv->dpos.at(idEpoch);
            recv->dpos.at(idEpoch) = Vector3d(p);
            const Double dr = (recv->dpos.at(idEpoch)-posOld).r();
            if(dr > maxChangePos)
            {
              maxChangePos     = dr;
              infoMaxChangePos = "  "+recv->name()+": "+recv->timeCorrected(idEpoch).dateTimeStr()+" posChange      = "+(1e3*dr)%"%6.1f"s+" mm";
            }
          }
        }

        // update wet troposphere
        // ----------------------
        if(recv->indexParameterTropoWet)
        {
          recv->xTropoWet += x.row(normalEquationInfo.index(recv->indexParameterTropoWet), parametrizationTropoWet->parameterCount());
          std::vector<UInt>   index;
          std::vector<Double> factor;
          for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
          {
            parametrizationTropoWet->factors(std::max(recv->timeCorrected(idEpoch), times.at(0)), index, factor);
            Double z = 0;
            for(UInt k=0; k<factor.size(); k++)
              z += factor.at(k) * recv->xTropoWet(index.at(k));
            const Double zOld = recv->zenitDelayWet.at(idEpoch);
            recv->zenitDelayWet.at(idEpoch) = z;
            if(std::fabs(z-zOld)> maxChangeWet)
            {
              maxChangeWet     = std::fabs(z-zOld);
              infoMaxChangeWet = "  "+recv->name()+": "+recv->timeCorrected(idEpoch).dateTimeStr()+" wetChange      = "+(1e3*(z-zOld))%"%6.1f"s+" mm";
            }
          }
        }

        // update troposphere gradient
        // ---------------------------
        if(recv->indexParameterTropoGradient)
        {
          recv->xTropoGradient += x.row(normalEquationInfo.index(recv->indexParameterTropoGradient), 2*parametrizationTropoGradient->parameterCount());
          std::vector<UInt>   index;
          std::vector<Double> factor;
          for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
          {
            parametrizationTropoGradient->factors(std::max(recv->timeCorrected(idEpoch), times.at(0)), index, factor);
            Double dx  = 0,  dy = 0;
            for(UInt k=0; k<factor.size(); k++)
            {
              dx += factor.at(k) * recv->xTropoGradient(2*index.at(k)+0);
              dy += factor.at(k) * recv->xTropoGradient(2*index.at(k)+1);
            }
            const Double dxOld = recv->gradientX.at(idEpoch);
            const Double dyOld = recv->gradientY.at(idEpoch);
            recv->gradientX.at(idEpoch) = dx;
            recv->gradientY.at(idEpoch) = dy;
            const Double dxdy = sqrt(pow(dx-dxOld,2) + pow(dy-dyOld,2));
            if(dxdy > maxChangeGradient)
            {
              maxChangeGradient     = dxdy;
              infoMaxChangeGradient = "  "+recv->name()+": "+recv->timeCorrected(idEpoch).dateTimeStr()+" gradientChange = "+(1e3*(dx-dxOld))%"%6.1f"s+", "+(1e3*(dy-dyOld))%"%.1f"s+" mm";
            }
          }
        }
      } // for(idRecv)

    // ===========================================================

    // update signal biases (at all nodes)
    // -----------------------------------
    for(auto &recv : receiver)
      if(recv->indexParameterBias)
      {
        const Vector dBias = recv->Bias * x.row(normalEquationInfo.index(recv->indexParameterBias), recv->Bias.columns());
        for(UInt idType=0; idType<dBias.size(); idType++)
          recv->signalBias.bias.at(idType) += dBias(idType);
        for(UInt idType=0; idType<dBias.size(); idType++)
          if(std::fabs(dBias(idType)) > maxChangeBias)
          {
            maxChangeBias = std::fabs(dBias(idType));
            infoMaxChangeBias = "  "+recv->name()+": "+recv->signalBias.type.at(idType).str()+"              biasChange     = "+(1e3*dBias(idType))%"%6.1f mm"s;
          }
      }

    // ===========================================================

    // compute new net translation/rotation
    // ------------------------------------
    if((normalEquationInfo.estimationType & Gnss::NormalEquationInfo::CONSTRAINT_RECEIVER_POSITION_NONETTRANSLATION) ||
       (normalEquationInfo.estimationType & Gnss::NormalEquationInfo::CONSTRAINT_RECEIVER_POSITION_NONETROTATION))
    {
      netRotation    = Vector(3);
      netTranslation = Vector(3);
      UInt count = 0;
      for(UInt idRecv=0; idRecv<receiver.size(); idRecv++)
        if(noNetStations(idRecv))
          for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
            if(receiver.at(idRecv)->useable(idEpoch))
            {
              const Vector   dpos = receiver.at(idRecv)->dpos.at(idEpoch).vector();
              const Vector3d pos  = receiver.at(idRecv)->approxPosition * (1/DEFAULT_R);
              Matrix R(3,3);
              R(0,1) =  pos.z(); R(0,2) = -pos.y();
              R(1,0) = -pos.z(); R(1,2) =  pos.x();
              R(2,0) =  pos.y(); R(2,1) = -pos.x();
              netRotation    += R * dpos;
              netTranslation += dpos;
              count++;
            }
      Parallel::reduceSum(count,          0, normalEquationInfo.comm);
      Parallel::reduceSum(netRotation,    0, normalEquationInfo.comm);
      Parallel::reduceSum(netTranslation, 0, normalEquationInfo.comm);
      if(Parallel::isMaster(normalEquationInfo.comm))
      {
        netRotation    *= 1./count;
        netTranslation *= 1./count;
        if(count == 0)
        {
          logWarning<<"No stations for no net translation/rotation computation"<<Log::endl;
          netRotation = netTranslation = Vector(3);
        }
      }
      Parallel::broadCast(netRotation,    0, normalEquationInfo.comm);
      Parallel::broadCast(netTranslation, 0, normalEquationInfo.comm);

      // info about net translation/rotation
      if(Parallel::isMaster(normalEquationInfo.comm) && (norm(netTranslation)>1e-5))
      {
        logInfo<<"                            netTranslation = "
              <<(1e3*netTranslation(0))%"%6.1f, "s
              <<(1e3*netTranslation(1))%"%6.1f, "s
              <<(1e3*netTranslation(2))%"%6.1f mm"s<<Log::endl;
      }

      if(Parallel::isMaster(normalEquationInfo.comm) && (norm(netRotation)>1e-5))
      {
        logInfo<<"                            netRotation    = "
              <<(1e3*netRotation(0))%"%6.1f, "s
              <<(1e3*netRotation(1))%"%6.1f, "s
              <<(1e3*netRotation(2))%"%6.1f mm"s<<Log::endl;
      }
    } // if(!SINGLERECEIVER)

    // ===========================================================

    Gnss::checkMaxChange(maxChangePos,      infoMaxChangePos,      printStatistics, normalEquationInfo.comm);
    Gnss::checkMaxChange(maxChangeClock,    infoMaxChangeClock,    printStatistics, normalEquationInfo.comm);
    Gnss::checkMaxChange(maxChangeWet,      infoMaxChangeWet,      printStatistics, normalEquationInfo.comm);
    Gnss::checkMaxChange(maxChangeGradient, infoMaxChangeGradient, printStatistics, normalEquationInfo.comm);
    Gnss::checkMaxChange(maxChangeBias,     infoMaxChangeBias,     printStatistics, normalEquationInfo.comm);

    return std::max(maxChangeClock, std::max(maxChangePos, maxChangeBias));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationReceiverStationNetwork::writeResults(const Gnss::NormalEquationInfo &normalEquationInfo, const std::string &suffix)
{
  try
  {
    VariableList fileNameVariableList;
    addVariable("station", fileNameVariableList);

    Vector useable(receiver.size());
    for(UInt idRecv=0; idRecv<receiver.size(); idRecv++)
      if(normalEquationInfo.estimateReceiver.at(receiver.at(idRecv)->idRecv()) && receiver.at(idRecv)->useable())
        useable(idRecv) = 1.;
    Parallel::reduceSum(useable, 0, normalEquationInfo.comm);

    // ===========================================================

    if(!fileNameUsedStationList.empty() && Parallel::isMaster(normalEquationInfo.comm))
    {
      logStatus<<"write used station list to file <"<<fileNameUsedStationList(fileNameVariableList).appendBaseName(suffix)<<">"<<Log::endl;
      std::vector<std::string> usedStationList;
      for(UInt idRecv=0; idRecv<receiver.size(); idRecv++)
        if(useable(idRecv))
          usedStationList.push_back(receiver.at(idRecv)->name());
      writeFileStringList(fileNameUsedStationList(fileNameVariableList).appendBaseName(suffix), usedStationList);
    }

    // ===========================================================

    if(!fileNameGrid.empty())
    {
      // collect positions
      Matrix data(receiver.size(), 6); // x, y, z , dx, dy, dz
      for(UInt idRecv=0; idRecv<receiver.size(); idRecv++)
        if(normalEquationInfo.estimateReceiver.at(receiver.at(idRecv)->idRecv()) && receiver.at(idRecv)->useable())
        {
          // Build average to get dpos
          const Vector3d dpos = std::accumulate(receiver.at(idRecv)->dpos.begin(), receiver.at(idRecv)->dpos.end(), Vector3d()) / receiver.at(idRecv)->dpos.size();
          // Get exact position
          const Vector3d pos = receiver.at(idRecv)->approxPosition + dpos;
          // dpos in NEU (north, east, up) system
          const Vector3d dposNEU = receiver.at(idRecv)->lnof2trf.inverseTransform(dpos);
          copy(pos.vector().trans(),     data.slice(idRecv, 0, 1, 3));
          copy(dposNEU.vector().trans(), data.slice(idRecv, 3, 1, 3));
        }
      Parallel::reduceSum(data, 0, normalEquationInfo.comm);

      if(Parallel::isMaster(normalEquationInfo.comm))
      {
        logStatus<<"write receiver grid file <"<<fileNameGrid(fileNameVariableList).appendBaseName(suffix)<<">"<<Log::endl;
        std::vector<Vector3d>            point;
        std::vector<std::vector<Double>> value(3);
        for(UInt idRecv=0; idRecv<receiver.size(); idRecv++)
          if(useable(idRecv))
          {
            point.push_back(Vector3d(data.slice(idRecv, 0, 1, 3)));
            value.at(0).push_back(data(idRecv, 3)); // north
            value.at(1).push_back(data(idRecv, 4)); // east
            value.at(2).push_back(data(idRecv, 5)); // up
          }
        writeFileGriddedData(fileNameGrid(fileNameVariableList).appendBaseName(suffix), GriddedData(Ellipsoid(), point, std::vector<Double>(point.size(), 1.), value));
      }
    }

    // ===========================================================

    if(!fileNamePosition.empty())
    {
      fileNameVariableList["station"]->setValue("****");
      logStatus<<"write receiver time series to files <"<<fileNamePosition(fileNameVariableList).appendBaseName(suffix)<<">"<<Log::endl;
      for(UInt idRecv=0; idRecv<receiver.size(); idRecv++)
        if(normalEquationInfo.estimateReceiver.at(receiver.at(idRecv)->idRecv()) && receiver.at(idRecv)->useable())
        {
          Vector3dArc   arc;
          Vector3dEpoch epoch;
          epoch.time     = times.at(static_cast<UInt>(std::round(0.5*times.size())));
          epoch.vector3d = receiver.at(idRecv)->approxPosition + std::accumulate(receiver.at(idRecv)->dpos.begin(), receiver.at(idRecv)->dpos.end(), Vector3d()) / receiver.at(idRecv)->dpos.size();
          arc.push_back(epoch);
          fileNameVariableList["station"]->setValue(receiver.at(idRecv)->stationName);
          InstrumentFile::write(fileNamePosition(fileNameVariableList).appendBaseName(suffix), arc);
        }
    }

    // ===========================================================

    if(!fileNamePositionSeries.empty())
    {
      fileNameVariableList["station"]->setValue("****");
      logStatus<<"write receiver time series to files <"<<fileNamePositionSeries(fileNameVariableList).appendBaseName(suffix)<<">"<<Log::endl;
      for(UInt idRecv=0; idRecv<receiver.size(); idRecv++)
        if(normalEquationInfo.estimateReceiver.at(receiver.at(idRecv)->idRecv()) && receiver.at(idRecv)->useable())
        {
          Matrix A(normalEquationInfo.idEpochs.size(), 11);
          UInt i = 0;
          for(UInt idEpoch : normalEquationInfo.idEpochs)
          {
            const Vector3d p1 = receiver.at(idRecv)->lnof2trf.inverseTransform(receiver.at(idRecv)->dpos.at(idEpoch));
            const Vector3d p2 = receiver.at(idRecv)->lnof2trf.inverseTransform(receiver.at(idRecv)->posDisplacement.at(idEpoch));
            A(i, 0) = receiver.at(idRecv)->timeCorrected(idEpoch).mjd();
            A(i, 1) = p1.x(); // north
            A(i, 2) = p1.y(); // east
            A(i, 3) = p1.z(); // up
            A(i, 4) = receiver.at(idRecv)->clockError(idEpoch)*LIGHT_VELOCITY;
            A(i, 5) = p2.x(); // north
            A(i, 6) = p2.y(); // east
            A(i, 7) = p2.z(); // up
            A(i, 8) = receiver.at(idRecv)->zenitDelayWet.at(idEpoch);
            A(i, 9) = receiver.at(idRecv)->gradientX.at(idEpoch);
            A(i,10) = receiver.at(idRecv)->gradientY.at(idEpoch);
            i++;
          }

          fileNameVariableList["station"]->setValue(receiver.at(idRecv)->stationName);
          writeFileMatrix(fileNamePositionSeries(fileNameVariableList).appendBaseName(suffix), A);
        } // for(idRecv)
    } // for(fileNamePositionSeries)

    // ===========================================================

    // write residuals
    // ---------------
    if(!fileNameResiduals.empty())
    {
      fileNameVariableList["station"]->setValue("****");
      logStatus<<"write receiver residuals to files <"<<fileNameResiduals(fileNameVariableList).appendBaseName(suffix)<<">"<<Log::endl;

      for(UInt idRecv=0; idRecv<receiver.size(); idRecv++)
        if(normalEquationInfo.estimateReceiver.at(receiver.at(idRecv)->idRecv()) && receiver.at(idRecv)->useable())
        {
          fileNameVariableList["station"]->setValue(receiver.at(idRecv)->stationName);
          InstrumentFile::write(fileNameResiduals(fileNameVariableList).appendBaseName(suffix), receiver.at(idRecv)->residuals(normalEquationInfo.idEpochs));
        } // for(idRecv)
    }

    // ===========================================================

    // Signal bias
    // -----------
    if(!fileNameOutSignal.empty())
    {
      fileNameVariableList["station"]->setValue("****");
      logStatus<<"write receiver signal bias to files <"<<fileNameOutSignal(fileNameVariableList).appendBaseName(suffix)<<">"<<Log::endl;
      for(auto &recv : receiver)
        if(normalEquationInfo.estimateReceiver.at(recv->idRecv()) && recv->useable())
        {
          GnssSignalBias signalBias = recv->signalBias;
          for(UInt idType=0; idType<signalBias.type.size(); idType++)
            if(signalBias.type.at(idType) == GnssType::PHASE)
              signalBias.bias.at(idType) = std::remainder(signalBias.bias.at(idType), LIGHT_VELOCITY/signalBias.type.at(idType).frequency());
          fileNameVariableList["station"]->setValue(recv->stationName);
          writeFileGnssSignalBias(fileNameOutSignal(fileNameVariableList).appendBaseName(suffix), signalBias);
        }
    }

    // ===========================================================

    // Estimated clocks
    // ----------------
    if(!fileNameOutClock.empty())
    {
      fileNameVariableList["station"]->setValue("****");
      logStatus<<"write receiver clock to files <"<<fileNameOutClock(fileNameVariableList).appendBaseName(suffix)<<">"<<Log::endl;

      for(auto &recv : receiver)
        if(normalEquationInfo.estimateReceiver.at(recv->idRecv()) && recv->useable())
        {
          fileNameVariableList["station"]->setValue(recv->stationName);

          MiscValueArc arc;
          for(UInt idEpoch : normalEquationInfo.idEpochs)
            if(recv->useable(idEpoch))
            {
              const Double clk = recv->clockError(idEpoch);
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

    // Troposphere values
    // ----------------
    if(!fileNameOutTropo.empty())
    {
      fileNameVariableList["station"]->setValue("****");
      logStatus<<"write troposphere values to files <"<<fileNameOutTropo(fileNameVariableList).appendBaseName(suffix)<<">"<<Log::endl;
      for(auto &recv : receiver)
        if(normalEquationInfo.estimateReceiver.at(recv->idRecv()) && recv->useable())
        {
          Matrix A(normalEquationInfo.idEpochs.size(), 12);
          for(UInt i=0; i<normalEquationInfo.idEpochs.size(); i++)
          {
            const UInt idEpoch = normalEquationInfo.idEpochs.at(i);
            Double zenithWetDelay, zenithDryDelay, gradientWetNorth, gradientDryNorth, gradientWetEast, gradientDryEast, aDry, aWet;
            troposphere->getAprioriValues(times.at(idEpoch), recv->idTropo, zenithDryDelay, zenithWetDelay,
                                          gradientDryNorth, gradientWetNorth, gradientDryEast, gradientWetEast, aDry, aWet);
            A(i, 0) = times.at(idEpoch).mjd();
            A(i, 1) = zenithDryDelay;                                     // tropospheric zenith dry delay [m] (only from model)
            A(i, 2) = zenithWetDelay   + recv->zenitDelayWet.at(idEpoch); // tropospheric zenith wet delay [m] (model + delta estimate)
            A(i, 3) = gradientDryNorth + recv->gradientX.at(idEpoch);     // tropospheric dry gradient - north direction [m] (model + delta estimate, due to same mapping function)
            A(i, 4) = gradientWetNorth;                                   // tropospheric wet gradient - north direction [m] (only from model)
            A(i, 5) = gradientDryEast  + recv->gradientY.at(idEpoch);     // tropospheric dry gradient - east component [m] (model + delta estimate, due to same mapping function)
            A(i, 6) = gradientWetEast;                                    // tropospheric wet gradient - east component [m] (only from model)
            A(i, 7) = recv->zenitDelayWet.at(idEpoch);                    // tropospheric zenith wet delay [m] (delta estimate)
            A(i, 8) = recv->gradientX.at(idEpoch);                        // tropospheric gradient - north [m] (delta estimate)
            A(i, 9) = recv->gradientY.at(idEpoch);                        // tropospheric gradient - east  [m] (delta estimate)
            A(i, 10) = aDry;                                              // dry mapping function coefficient a []
            A(i, 11) = aWet;                                              // wet mapping function coefficient a []
          }

          fileNameVariableList["station"]->setValue(recv->stationName);
          writeFileMatrix(fileNameOutTropo(fileNameVariableList).appendBaseName(suffix), A);
        }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

Bool GnssParametrizationReceiverStationNetwork::Receiver::isEpochEstimable(Gnss::AnalysisType analysisType, UInt idEpoch) const
{
  try
  {
    if(!useable(idEpoch))
      return FALSE;

    std::set<GnssType> systems;
    std::vector<GnssType> types;
    UInt count = 0;
    for(UInt idTrans=0; idTrans<gnss().transmitter.size(); idTrans++)
      if(gnss().transmitter.at(idTrans)->useable(idEpoch) && observation(idTrans, idEpoch) && observation(idTrans, idEpoch)->observationList(analysisType, types))
      {
        systems.insert(types.at(0) & GnssType::SYSTEM);
        count++;
      }

    return (count > (base->isKinematicPosition ? 3 : 0) + systems.size()); // position + 1 clock per system
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool GnssParametrizationReceiverStationNetwork::Receiver::isDesignMatrixReceiver(const Gnss::NormalEquationInfo &/*normalEquationInfo*/, UInt /*idTrans*/, UInt idEpoch) const
{
  try
  {
    if(!useable(idEpoch))
      return FALSE;

    if(indexParameterEpoch.at(idEpoch) ||
       indexParameterPosition          ||
       indexParameterTropoWet          ||
       indexParameterTropoGradient     ||
       indexParameterBias)
      return TRUE;

    return FALSE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationReceiverStationNetwork::Receiver::designMatrixReceiver(const Gnss::NormalEquationInfo &/*normalEquationInfo*/, const Gnss::ObservationEquation &eqn, Gnss::DesignMatrix &A) const
{
  try
  {
    // Epoch parameter
    // ---------------
    if(indexParameterEpoch.at(eqn.idEpoch))
    {
      MatrixSlice Design(A.column(indexParameterEpoch.at(eqn.idEpoch)));
      copy(eqn.A.column(Gnss::ObservationEquation::idxClockRecv,1), Design.column(0)); // clock error
      if(base->isKinematicPosition) // kinematic position
        matMult(1., eqn.A.column(Gnss::ObservationEquation::idxPosRecv,3), gnss().earthRotation->rotaryMatrix(eqn.timeRecv).matrix().trans(), Design.column(1, 3));
    } // if(indexParameterEpoch)

    // ===========================================================

    const Time timeRecv = std::max(eqn.timeRecv, base->times.at(0));

    // position
    // --------
    if((indexParameterPosition) && base->parametrizationPosition->parameterCount())
    {
      const Matrix B = eqn.A.column(Gnss::ObservationEquation::idxPosRecv,3) * gnss().earthRotation->rotaryMatrix(eqn.timeRecv).matrix().trans();
      designMatrixTemporal(base->parametrizationPosition, timeRecv, B, indexParameterPosition, A);
    }

    // troposphere wet
    // ---------------
    if(indexParameterTropoWet)
    {
      const Double mappingFunctionWet = base->troposphere->mappingFunctionWet(timeRecv, idTropo, eqn.azimutRecvLocal, eqn.elevationRecvLocal);
      const Matrix B = mappingFunctionWet * eqn.A.column(Gnss::ObservationEquation::idxRange,1);
      designMatrixTemporal(base->parametrizationTropoWet, timeRecv, B, indexParameterTropoWet, A);
    }

    // troposphere gradient
    // --------------------
    if(indexParameterTropoGradient)
    {
      Double dx, dy;
      base->troposphere->mappingFunctionGradient(timeRecv, idTropo, eqn.azimutRecvLocal, eqn.elevationRecvLocal, dx, dy);
      Matrix B(eqn.A.rows(), 2);
      axpy(dx, eqn.A.column(Gnss::ObservationEquation::idxRange,1), B.column(0));
      axpy(dy, eqn.A.column(Gnss::ObservationEquation::idxRange,1), B.column(1));
      designMatrixTemporal(base->parametrizationTropoGradient, timeRecv, B, indexParameterTropoGradient, A);
    }

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

void GnssParametrizationReceiverStationNetwork::Receiver::designMatrixTemporal(ParametrizationTemporalPtr parametrization, const Time &time, const_MatrixSliceRef B,
                                                                               const Gnss::ParameterIndex &index, Gnss::DesignMatrix &A)
{
  try
  {
    std::vector<UInt>   idx;
    std::vector<Double> factor;
    parametrization->factors(time, idx, factor);
    MatrixSlice Design(A.column(index));
    for(UInt i=0; i<factor.size(); i++)
      axpy(factor.at(i), B, Design.column(B.columns()*idx.at(i), B.columns()));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string GnssParametrizationReceiverStationNetwork::Receiver::name() const
{
  return stationName;
}

/***********************************************/

Vector3d GnssParametrizationReceiverStationNetwork::Receiver::position(UInt idEpoch) const
{
  const Rotary3d rotEarth = gnss().earthRotation->rotaryMatrix(timeCorrected(idEpoch));
  return rotEarth.inverseRotate( approxPosition + posDisplacement.at(idEpoch) + dpos.at(idEpoch) );
}

/***********************************************/

Vector3d GnssParametrizationReceiverStationNetwork::Receiver::velocity(UInt idEpoch) const
{
  return crossProduct(Vector3d(0.,0.,7.29211585531e-5), position(idEpoch));
}

/***********************************************/

Transform3d GnssParametrizationReceiverStationNetwork::Receiver::celestial2localFrame(UInt idEpoch) const
{
  const Rotary3d rotEarth = gnss().earthRotation->rotaryMatrix(timeCorrected(idEpoch));
  return inverse(lnof2trf) * rotEarth;
}

/***********************************************/

Transform3d GnssParametrizationReceiverStationNetwork::Receiver::local2antennaFrame(UInt idEpoch) const
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

Double GnssParametrizationReceiverStationNetwork::Receiver::troposphere(UInt idEpoch, Angle azimut, Angle elevation) const
{
  if(!base->troposphere)
    return 0;

  const Time t = std::max(timeCorrected(idEpoch), base->times.at(0));

  // apriori value
  Double delay = base->troposphere->slantDelay(t, idTropo, azimut, elevation);

  // estimated wet effect
  delay += base->troposphere->mappingFunctionWet(t, idTropo, azimut, elevation) * zenitDelayWet.at(idEpoch);

  // estimated gradient
  Double dx, dy;
  base->troposphere->mappingFunctionGradient(t, idTropo, azimut, elevation, dx, dy);
  delay += dx*gradientX.at(idEpoch) + dy*gradientY.at(idEpoch);

  return delay;
}

/***********************************************/

std::vector<GnssType> GnssParametrizationReceiverStationNetwork::Receiver::definedTypes(UInt idEpoch) const
{
  try
  {
    UInt idRecv = stationInfo.findReceiver(base->times.at(idEpoch));
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

Vector GnssParametrizationReceiverStationNetwork::Receiver::antennaVariations(UInt idEpoch, Angle azimut, Angle elevation, const std::vector<GnssType> &types) const
{
  try
  {
    Vector corr(types.size());
    corr += stationInfo.antennaVariations(timeCorrected(idEpoch), azimut, elevation, types, base->noPatternFoundAction);
    corr += signalBias.compute(types); // Code/Phase biases

    return corr;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector GnssParametrizationReceiverStationNetwork::Receiver::accuracy(UInt idEpoch, Angle azimut, Angle elevation, const std::vector<GnssType> &types) const
{
  try
  {
    return stationInfo.accuracy(timeCorrected(idEpoch), azimut, elevation, types, base->noPatternFoundAction);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Angle GnssParametrizationReceiverStationNetwork::Receiver::elevationCutOff() const
{
  return base->elevationCutOff;
}

/***********************************************/

Bool GnssParametrizationReceiverStationNetwork::Receiver::supportsIntegerAmbiguities(const Gnss::NormalEquationInfo &/*normalEquationInfo*/) const
{
  return base->integerAmbiguities;
}

/***********************************************/

Bool GnssParametrizationReceiverStationNetwork::Receiver::isCodeBiasEstimated(const Gnss::NormalEquationInfo &normalEquationInfo) const
{
  try
  {
    return base->estimateCodeBias && (normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_RECEIVER_SIGNALBIAS);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool GnssParametrizationReceiverStationNetwork::Receiver::isPhaseBiasEstimated(const Gnss::NormalEquationInfo &normalEquationInfo) const
{
  try
  {
    return base->estimatePhaseBias && (normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_RECEIVER_SIGNALBIAS);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
