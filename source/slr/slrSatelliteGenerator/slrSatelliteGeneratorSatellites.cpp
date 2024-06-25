/***********************************************/
/**
* @file slrSatelliteGeneratorSatellites.cpp
*
* @brief Provides a list of SLR satellites.
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#include "base/import.h"
#include "base/polynomial.h"
#include "inputOutput/logging.h"
#include "files/filePlatform.h"
#include "files/fileInstrument.h"
#include "files/fileStringTable.h"
#include "slr/slrSatellite.h"
#include "slr/slrSatelliteGenerator/slrSatelliteGeneratorSatellites.h"

/***********************************************/

SlrSatelliteGeneratorSatellites::SlrSatelliteGeneratorSatellites(Config &config)
{
  try
  {
    readConfig(config, "inputfileSatelliteList", fileNamesSatelliteList, Config::MUSTSET,  "{groopsDataDir}/slr/satellites/satelliteList.txt", "ascii file with satellite names, used to loop variable {satellite}");
    readConfig(config, "inputfileSatelliteInfo", fileNameSatelliteInfo,  Config::MUSTSET,  "{groopsDataDir}/slr/satellites/satelliteInfo/satelliteInfo.{satellite}.xml", "variable {satellite} available");
    readConfig(config, "inputfileOrbit",         fileNameOrbit,          Config::MUSTSET,  "{satellite}_orbit_{loopTime:%D}.dat",    "variable {satellite} available");
    readConfig(config, "inputfileAttitude",      fileNameAttitude,       Config::OPTIONAL, "{satellite}_attitude_{loopTime:%D}.dat", "variable {satellite} available");
    readConfig(config, "interpolationDegree",    interpolationDegree,    Config::DEFAULT,  "7", "for orbit interpolation and velocity calculation");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SlrSatelliteGeneratorSatellites::init(const std::vector<Time> &/*times*/, std::vector<SlrSatellitePtr> &satellites)
{
  try
  {
    logStatus<<"init satellites"<<Log::endl;
    std::vector<std::string> satelliteList;
    for(const auto &fileNameSatelliteList : fileNamesSatelliteList)
    {
      std::vector<std::string> list;
      readFileStringList(fileNameSatelliteList, list);
      satelliteList.insert(satelliteList.end(), list.begin(), list.end());
    }
    VariableList varList;
    for(const std::string &satellite : satelliteList)
    {
      try
      {
        varList.setVariable("satellite", satellite);

        Platform platform;
        readFilePlatform(fileNameSatelliteInfo(varList), platform);
        platform.name = satellite;

        // read orbit
        // ----------
        OrbitArc orbit;
        try
        {
          orbit = InstrumentFile::read(fileNameOrbit(varList));
        }
        catch(std::exception &/*e*/)
        {
          logWarningOnce<<"Unable to read orbit file <"<<fileNameOrbit(varList)<<">, disabling satellite."<<Log::endl;
          continue;
        }

        // read attitude
        // -------------
        StarCameraArc attitude;
        if(!fileNameAttitude.empty())
        {
          try
          {
            attitude = InstrumentFile::read(fileNameAttitude(varList));
          }
          catch(std::exception &/*e*/)
          {
            logWarningOnce<<"Unable to read attitude file <"<<fileNameAttitude(varList)<<">, disabling satellite."<<Log::endl;
            continue;
          }
        }
        if(!attitude.size()) // simulate along track, cross, nadir
        {
          for(UInt i=0; i<orbit.size(); i++)
          {
            StarCameraEpoch epoch;
            epoch.time   = orbit.at(i).time;
            epoch.rotary = Rotary3d(orbit.at(i).velocity, crossProduct(orbit.at(i).velocity, orbit.at(i).position));
            attitude.push_back(epoch);
          }
        }

        Arc::checkSynchronized({orbit, attitude});
        std::vector<Time> times = orbit.times();
        const Matrix data = orbit.matrix();
        Matrix pos     = data.column(1, 3);
        Matrix vel     = data.column(4, 3);
        Matrix srf2crf = attitude.matrix().column(1, 4);

        satellites.push_back(std::make_shared<SlrSatellite>(platform, times, pos, vel, srf2crf, interpolationDegree));
      }
      catch(std::exception &e)
      {
        logWarningOnce<<"Unable to initialize satellite: "<<e.what()<<Log::endl;
      }
    } // for(idSat)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
