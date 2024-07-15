/***********************************************/
/**
* @file sinexEccentricties2SlrPlatform.cpp
*
* @brief Platform file from SINEX eccentricties.
*
* @author Torsten Mayer-Guerr
* @date 2022-12-12
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Reads metadata like station name, station number, approximate station position and station eccentricities
from \href{https://ilrs.gsfc.nasa.gov/network/site_procedures/eccentricity.html}{Station Eccentricities Sinex File}
(une version) and write them to the \configFile{outputfileStationInfo}{platform} for each station.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "inputOutput/fileSinex.h"
#include "files/filePlatform.h"

/***** CLASS ***********************************/

/** @brief Platform file from SINEX eccentricties.
* @ingroup programsConversionGroup */
class SinexEccentricties2SlrPlatform
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(SinexEccentricties2SlrPlatform, SINGLEPROCESS, "Platform file from SINEX eccentricties.", Conversion, Slr)

/***********************************************/

void SinexEccentricties2SlrPlatform::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNamePlatform, fileNameSinex;
    std::vector<std::string> stationNames;
    std::string variableLoopStation;

    readConfig(config, "outputfileStationInfo", fileNamePlatform,    Config::MUSTSET,  "",        "loop variable is replaced with station name");
    readConfig(config, "variableLoopStation",   variableLoopStation, Config::DEFAULT,  "station", "variable name for station loop");
    readConfig(config, "inputfileSinex",        fileNameSinex,       Config::MUSTSET,  "",        "SINEX file (.snx or .ssc)");
    readConfig(config, "stationName",           stationNames,        Config::OPTIONAL, "",        "convert only these stations");
    if(isCreateSchema(config)) return;

    logStatus<<"read SINEX file <"<<fileNameSinex<<">"<<Log::endl;
    Sinex sinex;
    readFileSinex(fileNameSinex, sinex);

    // SITE/ECCENTRICITY
    // -----------------
    std::map<std::string, Platform> platforms;
    for(std::string &line : sinex.findBlock("SITE/ECCENTRICITY")->lines)
    {
      // *SITE PT SOLN T DATA_START__ DATA_END____ UNE UP______ NORTH___ EAST____        CDP-SOD_
      std::string name = String::lowerCase(String::trim(line.substr(1, 4)));
      if(stationNames.size() && std::find(stationNames.begin(), stationNames.end(), name) == stationNames.end())
        continue;
      // if(line.at(14) != 'L')
      //   continue;
      if(line.substr(42, 3) != "UNE")
        throw(Exception("SINEX file must provide eccentricities in UNE"));

      Platform &platform  = platforms[name];
      platform.markerName = name;
      auto station = std::make_shared<PlatformSlrStation>();
      station->name = name;
      if(line.size() >= 88)
        station->name = String::trim(line.substr(80, 8)); // CDP-SOD
      station->timeStart = Sinex::str2time(line, 16, FALSE);
      station->timeEnd   = Sinex::str2time(line, 29, TRUE);
      station->position = Vector3d(String::toDouble(line.substr(55, 8)), String::toDouble(line.substr(64, 8)), String::toDouble(line.substr(46, 8)));
      platform.equipments.push_back(station);
    }

    // SITE/ID
    // -------
    for(std::string &line : sinex.findBlock("SITE/ID")->lines)
    {
      // *         1         2         3         4         5         6         7         8
      // *12345678901234567890123456789012345678901234567890123456789012345678901234567890
      // *Code PT __DOMES__ T _STATION DESCRIPTION__ APPROX_LON_ APPROX_LAT_ _APP_H_     CDP-SOD_
      line.resize(std::max(line.size(), UInt(44)), ' ');
      std::string name = String::lowerCase(String::trim(line.substr(1, 4)));
      if(stationNames.size() && std::find(stationNames.begin(), stationNames.end(), name) == stationNames.end())
        continue;
      Platform &platform = platforms[name];
      platform.markerName   = name;
      platform.markerNumber = String::trim(line.substr(9, 9));
      platform.comment      = String::trim(line.substr(21, 11));
      // approxPosition
      if(line.size() >= 75)
      {
        const Double longitude   = String::toDouble(line.substr(44, 3)) + String::toDouble(line.substr(48, 2))/60 + String::toDouble(line.substr(51, 4))/3600;
        const Double latitude    = String::toDouble(line.substr(56, 3)) + (String::startsWith(String::trim(line.substr(56, 3)), "-") ? -1 : 1) * std::fabs(String::toDouble(line.substr(60, 2))/60 + String::toDouble(line.substr(62, 5))/3600);
        const Double height      = String::toDouble(line.substr(68, 7));
        platform.approxPosition = Ellipsoid()(Angle(longitude), Angle(latitude), height);
      }
      if(line.size() >= 88)
      {
        const std::string CDP_SOD = String::trim(line.substr(80, 8)); // CDP-SOD
        auto iter = std::find_if(platform.equipments.begin(), platform.equipments.end(),
                                 [&](const auto &x) {return (x->name == CDP_SOD);});
        if(iter != platform.equipments.end())
        {
          auto station = std::dynamic_pointer_cast<PlatformSlrStation>(*iter);
          station->serial  = String::trim(line.substr(9, 9));
          station->comment = String::trim(line.substr(32, 11));
        }
      }
    }

    // write platform files
    // --------------------
    VariableList varList;
    varList.setVariable(variableLoopStation, "****");
    logStatus<<"write platform files <"<<fileNamePlatform(varList)<<">"<<Log::endl;
    for(auto &platformPair : platforms)
    {
      Platform &platform = platformPair.second;
      if(!platform.equipments.size())
        logWarning<<"station "<<platform.markerName<<" without laser station"<<Log::endl;
      std::stable_sort(platform.equipments.begin(), platform.equipments.end(), [](auto &p1, auto &p2) {return p1->timeStart < p2->timeStart;});
      for(UInt i=1; i<platform.equipments.size(); i++)
        if(platform.equipments.at(i)->timeStart < platform.equipments.at(i-1)->timeEnd)
          logWarning<<"overlapping time intervals between "<<platform.equipments.at(i-1)->name<<" ("<<platform.equipments.at(i-1)->serial<<") "
                    <<"ENU("<<platform.equipments.at(i-1)->position.x()<<", "<<platform.equipments.at(i-1)->position.y()<<", "<<platform.equipments.at(i-1)->position.z()<<") ends "<<platform.equipments.at(i-1)->timeEnd.dateTimeStr()
                    <<" and "<<platform.equipments.at(i)->name<<" ("<<platform.equipments.at(i)->serial<<") "
                    <<"ENU("<<platform.equipments.at(i)->position.x()<<", "<<platform.equipments.at(i)->position.y()<<", "<<platform.equipments.at(i)->position.z()<<") starts "<<platform.equipments.at(i)->timeStart.dateTimeStr()<<Log::endl;
      varList.setVariable(variableLoopStation, platform.markerName);
      writeFilePlatform(fileNamePlatform(varList), platform);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
