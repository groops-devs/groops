/***********************************************/
/**
* @file viennaMappingFunctionStation2File.cpp
*
* @brief Converts VMF station time series to one file.
*
* @author Torsten Mayer-Guerr
* @date 2017-01-11
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Converts Vienna Mapping Functions (VMF) station time series into \file{GROOPS file format}{griddedDataTimeSeries}.

Station-wise VMF data for GNSS is available at: \url{https://vmf.geo.tuwien.ac.at/trop_products/GNSS/}
)";

/***********************************************/

#include "programs/program.h"
#include "base/griddedData.h"
#include "base/string.h"
#include "inputOutput/file.h"
#include "files/filePlatform.h"
#include "files/fileGriddedDataTimeSeries.h"
#include "files/fileGriddedData.h"

/***** CLASS ***********************************/

/** @brief Converts VMF station time series to one file.
* @ingroup programsConversionGroup */
class ViennaMappingFunctionStation2File
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(ViennaMappingFunctionStation2File, SINGLEPROCESS, "converts VMF station time series to one file.", Conversion)

/***********************************************/

void ViennaMappingFunctionStation2File::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameVmfCoefficients;
    FileName fileNameStationInfo, fileNameStation, fileNameIn;

    readConfig(config, "outputfileVmfCoefficients", fileNameVmfCoefficients, Config::OPTIONAL, "", "");
    readConfig(config, "inputfileStationInfo",      fileNameStationInfo,     Config::MUSTSET,  "{groopsDataDir}/gnss/receiverStation/stationInfo/igs/stationInfo.{station}.xml", "");
    readConfig(config, "inputfileStation",          fileNameStation,         Config::MUSTSET,  "", "");
    readConfig(config, "inputfileVmf",              fileNameIn,              Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    // ======================================================

    VariableList fileNameVariableList;
    addVariable("station", fileNameVariableList);

    logStatus<<"read station list"<<Log::endl;
    GriddedData grid;
    std::vector<std::string> stationNames;
    {
      Ellipsoid ellipsoid;
      InFile file(fileNameStation);
      std::string line;
      while(std::getline(file, line))
      {
        std::string name;
        Double      lat, lon, h;
        std::stringstream ss(line);
        ss>>name>>lat>>lon>>h;
        stationNames.push_back(String::lowerCase(name));
        grid.points.push_back(ellipsoid(Angle(lon*DEG2RAD), Angle(lat*DEG2RAD), h));
      }

      VariableList fileNameVariableList;
      addVariable("station", fileNameVariableList);
      for(UInt idStation=0; idStation<stationNames.size(); idStation++)
      {
        try
        {
          // read position from stationInfo
          Platform stationInfo;
          fileNameVariableList["station"]->setValue(stationNames.at(idStation));
          readFilePlatform(fileNameStationInfo(fileNameVariableList), stationInfo);

          if((stationInfo.approxPosition-grid.points.at(idStation)).r() > 1e3)
            logWarning<<stationNames.at(idStation)<<": r = "<<0.001*(stationInfo.approxPosition-grid.points.at(idStation)).r()<<" km"<<Log::endl;
          grid.points.at(idStation) = stationInfo.approxPosition;
        }
        catch(std::exception &)
        {
          logWarning<<stationNames.at(idStation)<<": stationInfo not found"<<Log::endl;
        }
      }
    }

    // ======================================================

    logStatus<<"read coefficients from files"<<Log::endl;
    std::vector<Time>   times;
    std::vector<Matrix> data; // ah, aw, zhd, zwd, gnh, geh, gnw, gew;
    InFile file(fileNameIn);
    for(;;)
    {
      std::string line;
      while(std::getline(file, line) && line.size() && (line.at(0) == '!' || line.at(0) == '#')); // skip header/comment lines
      if(line.empty())
        break;

      std::stringstream ss(line);
      std::string name;
      Double      mjd;
      Double      ah, aw, zhd, zwd;
      Double      pressure, temperature, waterVapour;
      Double      gnh = 0, geh = 0, gnw = 0, gew = 0;
      ss>>name>>mjd>>ah>>aw>>zhd>>zwd>>pressure>>temperature>>waterVapour>>gnh>>geh>>gnw>>gew;

      name = String::lowerCase(name);
      const UInt idStation = std::distance(stationNames.begin(), std::find(stationNames.begin(), stationNames.end(), name));
      if(idStation >= stationNames.size())
      {
        logWarning<<name<<" not in list -> use station info"<<Log::endl;
        Platform stationInfo; // read position from stationInfo
        try
        {
          fileNameVariableList["station"]->setValue(name);
          readFilePlatform(fileNameStationInfo(fileNameVariableList), stationInfo);
        }
        catch(std::exception &)
        {
          logWarning<<name<<": skip station without coordinates"<<Log::endl;
          continue;
        }

        stationNames.push_back(name);
        grid.points.push_back(stationInfo.approxPosition);
        for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
        {
          Matrix newData(stationNames.size(), 8, NAN_EXPR);
          std::swap(newData, data.at(idEpoch));
          copy(newData, data.at(idEpoch).row(0, newData.rows()));
        }
      }

      const UInt idEpoch = std::distance(times.begin(), std::find(times.begin(), times.end(), mjd2time(mjd)));
      if(idEpoch >= times.size())
      {
        times.push_back(mjd2time(mjd));
        data.resize(times.size(), Matrix(stationNames.size(), 8, NAN_EXPR));
      }

      data.at(idEpoch)(idStation, 0) = ah;
      data.at(idEpoch)(idStation, 1) = aw;
      data.at(idEpoch)(idStation, 2) = zhd;
      data.at(idEpoch)(idStation, 3) = zwd;
      data.at(idEpoch)(idStation, 4) = gnh/1000; // mm -> m
      data.at(idEpoch)(idStation, 5) = geh/1000;
      data.at(idEpoch)(idStation, 6) = gnw/1000;
      data.at(idEpoch)(idStation, 7) = gew/1000;
    }

    // ======================================================

    for(UInt idStation=0; idStation<stationNames.size(); idStation++)
    {
      UInt count = 0;
      for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
        if(!data.at(idEpoch)(idStation, 0))
          count++;
      if(count)
        logWarning<<stationNames.at(idStation)<<" missing "<<count<<" epochs of "<<times.size()<<Log::endl;
    }

    // ======================================================

    logStatus<<"write VMF coefficients to <"<<fileNameVmfCoefficients<<">"<<Log::endl;
    writeFileGriddedDataTimeSeries(fileNameVmfCoefficients, 1, times, grid, data);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/
