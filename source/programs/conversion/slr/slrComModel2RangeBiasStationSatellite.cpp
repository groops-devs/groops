/***********************************************/
/**
* @file slrComModel2RangeBiasStationSatellite.cpp
*
* @brief read CoM model.
*
* @author Torsten Mayer-Guerr
* @date 2022-12-27
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Converts the tables of CoM corrections of José Rodríguez
(\url{https://icts-yebes.oan.es/slr/com_models/models/}) into station/satellite
specific \configFile{outputfileRangeBias}{instrument}. Only the deviations to the default value in
\configFile{inputfileSatelliteInfo}{platform} are written. This program must be called for every
provided satellite. The range bias values can be used in
\configClass{parametrization:rangeBiasStationSatelliteApriori}{slrParametrizationType:rangeBiasStationSatelliteApriori}
in \program {SlrProcessing}.

Reference:
Rodriguez J., Otsubo T., Appleby G. Upgraded Modelling for the
Determination of Centre of Mass Corrections of Geodetic SLR
Satellites: Impact on Key Parameters of the Terrestrial Reference
Frame. Journal of Geodesy, 2019. doi: 10.1007/s00190-019-01315-0
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "inputOutput/file.h"
#include "files/filePlatform.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief read CoM model.
* @ingroup programsConversionGroup */
class SlrComModel2RangeBiasStationSatellite
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(SlrComModel2RangeBiasStationSatellite, SINGLEPROCESS, "read CoM model", Conversion, Slr)

/***********************************************/

void SlrComModel2RangeBiasStationSatellite::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName                 fileNameOut;
    FileName                 fileNamePlatform, fileNameIn;
    std::string              variableLoopStation;
    std::vector<std::string> stationNames;

    readConfig(config, "outputfileRangeBias",    fileNameOut,         Config::MUSTSET, "rangeBias.{station}.{satellite}.txt",  "MISCVALUE, variable {station} available");
    readConfig(config, "inputfileSatelliteInfo", fileNamePlatform,    Config::MUSTSET, "{groopsDataDir}/slr/satellites/satelliteInfo/satelliteInfo.{satellite}.xml", "");
    readConfig(config, "inputfile",              fileNameIn,          Config::MUSTSET, "com/com_<aji|et1|la2|las|lg1|lg2|str>.dat", "from Rodriguez model");
    readConfig(config, "variableLoopStation",    variableLoopStation, Config::DEFAULT,  "station", "variable name for station loop");
    readConfig(config, "stationName",            stationNames,        Config::OPTIONAL, "",        "convert only these stations");
    if(isCreateSchema(config)) return;

    logStatus<<"read platform file <"<<fileNamePlatform<<">"<<Log::endl;
    Platform platform;
    readFilePlatform(fileNamePlatform, platform);

    // read data
    // ---------
    std::map<std::string, MiscValueArc> arcs;     // per station
    std::map<std::string, Time>         timesEnd; // per station

    logStatus<<"read file <"<<fileNameIn<<">"<<Log::endl;
    InFile file(fileNameIn);
    std::string line;
    while(std::getline(file, line))
    {
      if(line.empty() || String::startsWith(line, "*"))
        continue;

      std::string station;
      UInt        dayStart, monthStart, yearStart;
      UInt        dayEnd, monthEnd, yearEnd;
      Double      wavelength, range;
      std::stringstream ss(line);
      ss>>station>>dayStart>>monthStart>>yearStart>>dayEnd>>monthEnd>>yearEnd>>wavelength>>range;

      if(stationNames.size() && std::find(stationNames.begin(), stationNames.end(), station) == stationNames.end())
        continue;

      const Time timeStart = date2time(yearStart, monthStart, dayStart);
      const Time timeEnd   = date2time(yearEnd,   monthEnd,   dayEnd);

      auto reflector = platform.findEquipment<PlatformLaserRetroReflector>(timeStart);
      if(!reflector)
       reflector = platform.findEquipment<PlatformLaserRetroReflector>(timeEnd);
      if(!reflector)
      {
        logWarning<<"reflector not found in ["<<timeStart.dateStr()<<", "<<timeEnd.dateStr()<<"] -> skipping"<<Log::endl;
        continue;
      }

      MiscValueEpoch epoch;
      epoch.time  = timeStart;
      epoch.value = 1e-3 * range - reflector->range(0, 0); // mm -> m
      arcs[station].push_back(epoch);
      timesEnd[station] = timeEnd;
    }

    // write results
    // -------------
    VariableList varList;
    varList.setVariable(variableLoopStation, "****");
    logStatus<<"write range biases to file <"<<fileNameOut(varList)<<">"<<Log::endl;
    for(auto &stationArcs : arcs)
    {
      varList.setVariable(variableLoopStation, stationArcs.first);
      if(stationArcs.second.size() && std::any_of(stationArcs.second.begin(), stationArcs.second.end(), [&](const auto &e){return (std::fabs(e.value) > 0.0001);}))
      {
        if(timesEnd[stationArcs.first] < date2time(2050, 1, 1))
        {
          MiscValueEpoch epoch;
          epoch.time  = timesEnd[stationArcs.first];
          epoch.value = 0;
          stationArcs.second.push_back(epoch);
        }
        InstrumentFile::write(fileNameOut(varList), stationArcs.second);
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
