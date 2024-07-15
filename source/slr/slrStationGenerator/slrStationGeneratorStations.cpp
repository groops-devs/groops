/***********************************************/
/**
* @file slrStationGeneratorStations.cpp
*
* @brief SLR ground station network.
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#include "base/import.h"
#include "base/string.h"
#include "base/planets.h"
#include "config/config.h"
#include "inputOutput/logging.h"
#include "inputOutput/system.h"
#include "files/filePlatform.h"
#include "files/fileInstrument.h"
#include "files/fileMatrix.h"
#include "files/fileStringTable.h"
#include "classes/gravityfield/gravityfield.h"
#include "classes/tides/tides.h"
#include "slr/slr.h"
#include "slr/slrStation.h"
#include "slr/slrSatellite.h"
#include "slr/slrStationGenerator/slrStationGenerator.h"
#include "slr/slrStationGenerator/slrStationGeneratorStations.h"

/***********************************************/

SlrStationGeneratorStations::SlrStationGeneratorStations(Config &config)
{
  try
  {
    readConfig(config, "inputfileStationList",               fileNameStationList,     Config::MUSTSET,  "{groopsDataDir}/slr/stations/stationList.txt", "ascii file with station names");
    readConfig(config, "inputfileStationInfo",               fileNameStationInfo,     Config::MUSTSET,  "{groopsDataDir}/slr/stations/stationInfo/stationInfo.{station}.xml", "station metadata");
    readConfig(config, "inputfileStationPosition",           fileNameStationPosition, Config::OPTIONAL, "{groopsDataDir}/slr/stations/position/slrf2020_CM/stationPosition.{station}.dat", "station position");
    readConfig(config, "inputfileObservations",              fileNameObs,             Config::OPTIONAL, "normalPoints_{satellite}_{station}_{loopTime:%D}.dat", "variable {station} {satellite} available");
    readConfig(config, "accuracy",                           accuracyExpr,            Config::MUSTSET,  "if(abs(residual)>30, NAN, accuracy)", "[m] used for weighting, variables: {residual}, {accuracy}, {redundancy}, {laserWavelength}, {azimut}, {elevation}");
    readConfig(config, "loadingDisplacement",                gravityfield,            Config::DEFAULT,  "",    "loading deformation");
    readConfig(config, "tidalDisplacement",                  tides,                   Config::DEFAULT,  "",    "tidal deformation");
    readConfig(config, "ephemerides",                        ephemerides,             Config::OPTIONAL, "jpl", "for tidal deformation");
    readConfig(config, "inputfileDeformationLoadLoveNumber", deformationName,         Config::MUSTSET,  "{groopsDataDir}/loading/deformationLoveNumbers_CM_Gegout97.txt", "");
    readConfig(config, "inputfilePotentialLoadLoveNumber",   potentialName,           Config::OPTIONAL, "{groopsDataDir}/loading/loadLoveNumbers_Gegout97.txt", "if full potential is given and not only loading potential");
    readConfig(config, "elevationCutOff",                    elevationCutOff,         Config::DEFAULT,  "5", "[degree] ignore observations below cutoff");
    readConfig(config, "interpolationDegree",                interpolationDegree,     Config::DEFAULT,  "7", "for position interpolation");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SlrStationGeneratorStations::init(const std::vector<Time> &times, const std::vector<SlrSatellitePtr> &satellites,
                                       EarthRotationPtr earthRotation, std::vector<SlrStationPtr> &stationsAll)
{
  try
  {
    logStatus<<"init station network"<<Log::endl;
    std::vector<std::string> stationNames;
    readFileStringList(fileNameStationList, stationNames);
    VariableList varList;
    std::vector<SlrStationPtr> stations;
    for(const std::string &stationName : stationNames)
    {
      try
      {
        varList.setVariable("station", stationName);

        Platform platform;
        readFilePlatform(fileNameStationInfo(varList), platform);
        platform.name = stationName;
        // is station available in this time interval?
        if(!std::any_of(times.begin(), times.end(), [&](const Time &time) {return platform.findEquipment<PlatformSlrStation>(time);}))
          continue;

        // approximate station position
        if(!fileNameStationPosition.empty())
        {
          try
          {
            Vector3dArc arc = InstrumentFile::read(fileNameStationPosition(varList));
            auto iter = (arc.size() == 1) ? arc.begin() : std::find_if(arc.begin(), arc.end(), [&](const Epoch &e){return e.time.isInInterval(times.front(), times.back());});
            if(iter != arc.end())
            {
              platform.approxPosition = iter->vector3d;
              // logInfo<<stationName<<" precise pos at "<<iter->time.dateTimeStr()<<Log::endl;
            }
            else
              throw(Exception(""));
          }
          catch(std::exception &/*e*/)
          {
            logInfo<<stationName<<" no precise pos !!!!!!!!!!!!"<<Log::endl;
          }
        }

        Matrix offset(times.size(), 3);
        for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
        {
          auto station = platform.findEquipment<PlatformSlrStation>(times.at(idEpoch));
          if(station)
            axpy(1., (station->position - platform.referencePoint(times.at(idEpoch))).vector().trans(), offset.row(idEpoch));
        }

        auto station = std::make_shared<SlrStation>(platform, times, platform.approxPosition, offset, accuracyExpr, interpolationDegree);

        // read observations
        // -----------------
        if(!fileNameObs.empty())
        {
          auto rotationCrf2Trf = std::bind(&EarthRotation::rotaryMatrix, earthRotation, std::placeholders::_1);
          station->observations.resize(satellites.size());
          for(UInt idSat=0; idSat<satellites.size(); idSat++)
          {
            VariableList varList;
            varList.setVariable("station",   stationName);
            varList.setVariable("satellite", satellites.at(idSat)->name());
            if(!System::exists(fileNameObs(varList)))
              continue;

            InstrumentFile file(fileNameObs(varList));
            for(UInt idPass=0; idPass<file.arcCount(); idPass++)
            {
              SatelliteLaserRangingArc arc = file.readArc(idPass);
              SlrObservationPtr obs = std::make_shared<SlrObservation>();
              obs->timesTrans       = arc.times();
              obs->observations     = arc.matrix().column(1);
              obs->sigmas0          = arc.matrix().column(2);
              obs->redundancies     = arc.matrix().column(3);
              obs->laserWavelength  = arc.matrix().column(5);
              if(obs->init(*station, *satellites.at(idSat), rotationCrf2Trf, elevationCutOff))
                station->observations.at(idSat).push_back(obs);
            }
          }

          if(std::none_of(station->observations.begin(), station->observations.end(), [](auto &o) {return o.size();}))
            continue; // station without observations
        }

        stations.push_back(station);
      }
      catch(std::exception &e)
      {
        logWarningOnce<<stationName<<" disabled: "<<e.what()<<Log::endl;
      }
    } // for(stationNames)

    // store valid stations
    // ---------------------
    stationsAll.insert(stationsAll.end(), stations.begin(), stations.end());
    logInfo<<"  "<<stations.size()<<" of "<<stationNames.size()<<" stations used"<<Log::endl;

    // tides & loading
    // ---------------
    logStatus<<"compute tides & loading"<<Log::endl;
    Vector hn, ln;
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

    std::vector<Vector3d> positions;
    for(auto &stat : stations)
      positions.push_back(stat->position(times.front()));

    Vector gravity(positions.size()); // normal gravity
    for(UInt i=0; i<gravity.size(); i++)
      gravity(i) = Planets::normalGravity(positions.at(i));

    std::vector<Rotary3d> rotEarth(times.size());
    for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
      rotEarth.at(idEpoch) = earthRotation->rotaryMatrix(times.at(idEpoch));

    std::vector<std::vector<Vector3d>> disp(positions.size(), std::vector<Vector3d>(times.size()));
    tides->deformation(times, positions, rotEarth, earthRotation, ephemerides, gravity, hn, ln, disp);
    gravityfield->deformation(times, positions, gravity, hn, ln, disp);
    tides        = nullptr;
    gravityfield = nullptr;

    // add displacements
    for(UInt i=0; i<stations.size(); i++)
      for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
        axpy(1., stations.at(i)->global2localFrame(times.at(idEpoch)).transform(disp.at(i).at(idEpoch)).vector().trans(), stations.at(i)->offset.row(idEpoch));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
