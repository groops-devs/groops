/***********************************************/
/**
* @file slrSinexDataHandling2Files.cpp
*
* @brief Convert SLR SINEX data handling file.
*
* @author Torsten Mayer-Guerr
* @date 2022-12-15
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Converts SLR range and time bias from \verb|ILRS_Data_Handling_File_xxxx.xx.xx.snx| provided at
\url{https://cddis.nasa.gov/archive/slr/products/resource/}. The range and time bias values can be used in
\configClass{parametrization:rangeBiasXxxApriori}{slrParametrizationType:rangeBiasStationApriori}
in \program {SlrProcessing}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "inputOutput/fileSinex.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Convert SLR SINEX data handling file.
* @ingroup programsConversionGroup */
class SlrSinexDataHandling2Files
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(SlrSinexDataHandling2Files, SINGLEPROCESS, "Convert SLR SINEX data handling file", Conversion, Slr)

/***********************************************/

void SlrSinexDataHandling2Files::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName                 fileNameRangeBiasStation, fileNameRangeBiasSatellite, fileNameTimeBias;
    FileName                 fileNameSinex, fileNameSatelliteId;
    std::string              variableLoopStation, variableLoopSatellite;
    std::vector<std::string> stationNames;

    readConfig(config, "outputfileRangeBiasStation",          fileNameRangeBiasStation,   Config::OPTIONAL, "rangeBias.{station}.txt",             "MISCVALUE [m]");
    readConfig(config, "outputfileRangeBiasStationSatellite", fileNameRangeBiasSatellite, Config::OPTIONAL, "rangeBias.{station}.{satellite}.txt", "MISCVALUE [m]");
    readConfig(config, "outputfileTimeBias",                  fileNameTimeBias,           Config::OPTIONAL, "timeBias.{station}.txt",              "MISCVALUES(bias [s], drift [s/d])");
    readConfig(config, "variableLoopStation",                 variableLoopStation,        Config::DEFAULT,  "station",                             "variable name for station loop");
    readConfig(config, "variableLoopSatellite",               variableLoopSatellite,      Config::DEFAULT,  "satellite",                           "variable name for satellite loop");
    readConfig(config, "inputfileSinex",                      fileNameSinex,              Config::MUSTSET,  "ILRS_Data_Handling_File.snx",         "SINEX file (.snx)");
    readConfig(config, "inputfileSatelliteId",                fileNameSatelliteId,        Config::OPTIONAL, "",                                    "table SP3 and satellite name");
    readConfig(config, "stationName",                         stationNames,               Config::OPTIONAL, "",                                    "convert only these stations");
    if(isCreateSchema(config)) return;

    std::map<std::string, MiscValueArc>                         rangeBiasesStation;
    std::map<std::string, std::map<std::string, MiscValueArc>>  rangeBiasesStationSatellite;
    std::map<std::string, MiscValuesArc>                        timeBiases;

    logStatus<<"read SINEX file <"<<fileNameSinex<<">"<<Log::endl;
    Sinex sinex;
    readFileSinex(fileNameSinex, sinex);

    // MODEL/RANGE_BIAS
    // ----------------------
    for(std::string &line : sinex.findBlock("MODEL/RANGE_BIAS")->lines)
    {
      // *         1         2         3         4         5         6         7         8
      // *12345678901234567890123456789012345678901234567890123456789012345678901234567890
      // *CODE PT UNIT T _DATA_START_ __DATA_END__ M __E-VALUE___ STD_DEV _E-RATE__ CMNTS
      line.resize(80, ' ');
      std::string station = String::lowerCase(String::trim(line.substr(1, 4)));
      if(stationNames.size() && std::find(stationNames.begin(), stationNames.end(), station) == stationNames.end())
        continue;
      if(line.at(42) != 'R') // range bias
        continue;

      Time timeStart = Sinex::str2time(line, 16, FALSE);
      Time timeEnd   = Sinex::str2time(line, 29, TRUE);
      MiscValueEpoch epoch;
      epoch.time  = timeStart;
      epoch.value = 1e-3 * String::toDouble(line.substr(44, 12)); // mm -> m

      std::string satId = String::trim(line.substr(6, 2));
      if(satId == "--")
        satId = "";

      if(satId.empty()) // station specific range bias
      {
        if(rangeBiasesStation[station].size() && (rangeBiasesStation[station].back().time >= timeStart-seconds2time(1)))
          rangeBiasesStation[station].remove(rangeBiasesStation[station].size()-1);
        rangeBiasesStation[station].push_back(epoch);
        // zero bias after end of interval
        epoch.time  = timeEnd;
        epoch.value = 0.;
        rangeBiasesStation[station].push_back(epoch);
      }
      else
      {
        if(rangeBiasesStationSatellite[station][satId].size() && (rangeBiasesStationSatellite[station][satId].back().time >= timeStart-seconds2time(1)))
          rangeBiasesStationSatellite[station][satId].remove(rangeBiasesStationSatellite[station][satId].size()-1);
        rangeBiasesStationSatellite[station][satId].push_back(epoch);
        // zero bias after end of interval
        epoch.time  = timeEnd;
        epoch.value = 0.;
        rangeBiasesStationSatellite[station][satId].push_back(epoch);
      }
    }

    // MODEL/TIME_BIAS
    // ---------------
    for(std::string &line : sinex.findBlock("MODEL/TIME_BIAS")->lines)
    {
      // *         1         2         3         4         5         6         7         8
      // *12345678901234567890123456789012345678901234567890123456789012345678901234567890
      // *CODE PT SOLN T START_DATE__ END_DATE____ M __E-VALUE___ STD_DEV _E-RATE__ UNIT   CMNTS
      line.resize(80, ' ');
      std::string station = String::lowerCase(String::trim(line.substr(1, 4)));
      if(stationNames.size() && std::find(stationNames.begin(), stationNames.end(), station) == stationNames.end())
        continue;
      if(line.at(42) != 'T')  // time bias
        continue;
      if((line.substr(6, 2) != "  ") && (line.substr(6, 2) != "--"))
      {
        logWarning<<"station-satellite specific time bias ignored in line '"<<line<<"'"<<Log::endl;
        continue;
      }

      Time timeStart = Sinex::str2time(line, 16, FALSE);
      Time timeEnd   = Sinex::str2time(line, 29, TRUE);
      if(timeBiases[station].size() && (timeBiases[station].back().time >= timeStart-seconds2time(1)))
        timeBiases[station].remove(timeBiases[station].size()-1);

      Double unit = 1.;
      if(String::trim(line.substr(75, 4)) == "ms")
        unit = 1e-3;
      else if(String::trim(line.substr(75, 4)) == "us")
        unit = 1e-6;
      else
      {
        logWarning<<"Unknown time unit (ms or us) '"<<line<<"'"<<Log::endl;
        continue;
      }

      MiscValuesEpoch epoch(2);
      epoch.time  = timeStart;
      epoch.values(0) = unit * String::toDouble(line.substr(44, 12)); // bias  ms/us -> s
      epoch.values(1) = unit * String::toDouble(line.substr(65, 9));  // drift ms/us -> s/d
      epoch.values(0) -= epoch.values(1) * 0.5 * (timeEnd-timeStart).mjd(); // move bias from mid to start of interval
      timeBiases[station].push_back(epoch);
      // zero bias after end of interval
      epoch.time  = timeEnd;
      epoch.values.setNull();
      timeBiases[station].push_back(epoch);
    }

    // write results
    // -------------
    if(!fileNameRangeBiasStation.empty())
    {
      VariableList varList;
      varList.setVariable(variableLoopStation, "****");
      logStatus<<"write station range biases to file <"<<fileNameRangeBiasStation(varList)<<">"<<Log::endl;
      for(auto &rangeBiasPair : rangeBiasesStation)
        if(rangeBiasPair.second.size())
        {
          varList.setVariable(variableLoopStation, rangeBiasPair.first);
          InstrumentFile::write(fileNameRangeBiasStation(varList), rangeBiasPair.second);
        }
    }

    if(!fileNameRangeBiasSatellite.empty())
    {
      VariableList varList;
      varList.setVariable(variableLoopStation,   "****");
      varList.setVariable(variableLoopSatellite, "****");
      logStatus<<"write station-satellite range biases to file <"<<fileNameRangeBiasSatellite(varList)<<">"<<Log::endl;
      for(auto &pairStation : rangeBiasesStationSatellite)
        for(auto &pairSatellite : pairStation.second)
          if(pairSatellite.second.size())
          {
            varList.setVariable(variableLoopStation,   pairStation.first);
            varList.setVariable(variableLoopSatellite, pairSatellite.first);
            InstrumentFile::write(fileNameRangeBiasSatellite(varList), pairSatellite.second);
          }
    }

    if(!fileNameTimeBias.empty())
    {
      VariableList varList;
      varList.setVariable(variableLoopStation, "****");
      logStatus<<"write time biases to file <"<<fileNameTimeBias(varList)<<">"<<Log::endl;
      for(auto &timeBiasPair : timeBiases)
        if(timeBiasPair.second.size())
        {
          varList.setVariable(variableLoopStation, timeBiasPair.first);
          InstrumentFile::write(fileNameTimeBias(varList), timeBiasPair.second);
        }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
