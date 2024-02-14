/***********************************************/
/**
* @file sinex2StationDiscontinuities.cpp
*
* @brief Convert station discontinuities from SINEX (e.g. ITRF20) to InstrumentMiscValue.
*
* A value of 1 means position discontinuity, a value of 2 means velocity discontinuity.
*
* @author Sebastian Strasser
* @date 2017-05-22
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Convert station discontinuities from
\href{http://www.iers.org/IERS/EN/Organization/AnalysisCoordinator/SinexFormat/sinex.html}{SINEX format}
(e.g. ITRF20) to \configFile{outputfileInstrument}{instrument} (MISCVALUE).
A value of 1 means position discontinuity, a value of 2 means velocity discontinuity.
Start and end epochs with value 0 are added in addition to the discontinuities from
SINEX to define continuity interval borders.

See also \program{Sinex2StationPosition} and \program{Sinex2StationPostSeismicDeformation}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "inputOutput/fileSinex.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Convert station discontinuities from SINEX (e.g. ITRF20) to InstrumentMiscValue.
* @ingroup programsConversionGroup */
class Sinex2StationDiscontinuities
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Sinex2StationDiscontinuities, SINGLEPROCESS, "Convert station discontinuities from SINEX (e.g. ITRF20) to InstrumentMiscValue.", Conversion, Gnss, Instrument)

/***********************************************/

void Sinex2StationDiscontinuities::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName outName, fileNameSinex;
    std::vector<std::string> stationNames;
    std::string variableLoopStation;

    readConfig(config, "outputfileInstrument",     outName,             Config::MUSTSET,  "",        "loop variable is replaced with station name (e.g. wtzz)");
    readConfig(config, "inputfileDiscontinuities", fileNameSinex,       Config::MUSTSET,  "",        "SINEX (e.g. ITRF20) station discontinuities");
    readConfig(config, "variableLoopStation",      variableLoopStation, Config::DEFAULT,  "station", "variable name for station loop");
    readConfig(config, "stationName",              stationNames,        Config::OPTIONAL, "",        "only export these stations");
    if(isCreateSchema(config)) return;

    logStatus<<"read input file <"<<fileNameSinex<<">"<<Log::endl;
    Sinex sinex;
    readFileSinex(fileNameSinex, sinex);
    const std::vector<std::string> &lines = sinex.findBlock("SOLUTION/DISCONTINUITY")->lines;
    std::map<std::string, MiscValueArc> stations;
    for(const auto &line : lines)
    {
      std::string name = String::lowerCase(String::trim(line.substr(1, 4)));
      if(stationNames.size() && std::find(stationNames.begin(), stationNames.end(), name) == stationNames.end())
        continue;

      Time time = Sinex::str2time(line, 16, TRUE);
      MiscValueEpoch epoch;
      epoch.time = time;
      std::string type = line.substr(42,1);
      if(type == "P")
        epoch.value = 1;
      else if(type == "V")
        epoch.value = 2;
      else
      {
        logWarning<<"Unknown discontinuity type for line: "<<line<<Log::endl;
        continue;
      }
      stations[name].push_back(epoch);
    }

    logStatus<<"write output files <"<<outName<<">"<<Log::endl;
    VariableList fileNameVariableList;
    for(auto &&station : stations)
    {
      station.second.sort();
      while(station.second.size() && station.second.back().time == date2time(2500, 1, 1))
        station.second.remove(station.second.size()-1);
      MiscValueEpoch epoch;
      station.second.insert(0, epoch); // start epoch
      epoch.time = date2time(2500, 1, 1);
      station.second.push_back(epoch); // end epoch

      fileNameVariableList.setVariable(variableLoopStation, station.first);
      InstrumentFile::write(outName(fileNameVariableList), station.second);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
