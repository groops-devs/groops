/***********************************************/
/**
* @file sinex2StationPosition.cpp
*
* @brief DEPRECATED. Please use Sinex2StationPositions instead.
*
* @author Sebastian Strasser
* @date 2016-12-07
*
* @deprecated Please use Sinex2StationPositions instead.
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
DEPRECATED. Please use \program{Sinex2StationPositions} instead.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "classes/timeSeries/timeSeries.h"
#include "inputOutput/fileSinex.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief DEPRECATED. Please use Sinex2StationPositions instead.
* @ingroup programsConversionGroup */
class Sinex2StationPosition
{
public:
  class Interval
  {
  public:
    std::string solutionId;
    Time        referenceTime, timeStart, timeEnd;
    Vector3d    position, velocity;

    Bool operator<(const Interval &other) const {return timeStart < other.timeStart;}
  };

  class Station
  {
  public:
    std::vector<Interval> intervals;
    std::vector<Time>     discontinuities;

    Vector3dArc arc(const std::vector<Time> &times, Bool extrapolateForward, Bool extrapolateBackward);
  };

  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Sinex2StationPosition, SINGLEPROCESS, "DEPRECATED. Please use Sinex2StationPositions instead.", Deprecated)
GROOPS_RENAMED_PROGRAM(GnssSinex2StationPosition, Sinex2StationPosition, date2time(2018, 12, 4))

/***********************************************/

void Sinex2StationPosition::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameInstrument, fileNameDiscontinuities, fileNameSinex;
    std::vector<std::string> stationNames;
    std::string variableLoopStation;
    TimeSeriesPtr timeSeries;
    Bool extrapolateForward, extrapolateBackward;

    readConfig(config, "outputfileInstrument",     fileNameInstrument,      Config::MUSTSET,   "",        "loop variable is replaced with station name (e.g. wtzz)");
    readConfig(config, "inputfileSinex",           fileNameSinex,           Config::MUSTSET,   "",        "SINEX file (.snx or .ssc)");
    readConfig(config, "inputfileDiscontinuities", fileNameDiscontinuities, Config::OPTIONAL,  "",        "discontinuities file per station; loop variable is replaced with station name (e.g. wtzz)");
    readConfig(config, "variableLoopStation",      variableLoopStation,     Config::DEFAULT,   "station", "variable name for station loop");
    readConfig(config, "stationName",              stationNames,            Config::OPTIONAL,  "",        "convert only these stations");
    readConfig(config, "timeSeries",               timeSeries,              Config::DEFAULT,   "",        "compute positions for these epochs based on velocity");
    readConfig(config, "extrapolateForward",       extrapolateForward,      Config::DEFAULT,   "0",       "also compute positions for epochs after last interval defined in SINEX file");
    readConfig(config, "extrapolateBackward",      extrapolateBackward,     Config::DEFAULT,   "0",       "also compute positions for epochs before first interval defined in SINEX file");
    if(isCreateSchema(config)) return;

    logWarning<<"DEPRECATED. Please use Sinex2StationPositions instead."<<Log::endl;

    logStatus<<"read SINEX file <"<<fileNameSinex<<">"<<Log::endl;
    std::map<std::string, Station> stations;
    Sinex sinex;
    readFileSinex(fileNameSinex, sinex);

    // SOLUTION/EPOCHS
    // ---------------
    for(std::string &line : sinex.findBlock("SOLUTION/EPOCHS")->lines)
    {
      std::string name = String::lowerCase(String::trim(line.substr(1, 4)));
      if(stationNames.size() && std::find(stationNames.begin(), stationNames.end(), name) == stationNames.end())
        continue;

      Interval interval;
      interval.solutionId = String::trim(line.substr(9, 4));
      interval.timeStart  = Sinex::str2time(line, 16, FALSE);
      interval.timeEnd    = Sinex::str2time(line, 29, TRUE);
      stations[name].intervals.push_back(interval);
    }
    for(auto &station : stations)
      std::sort(station.second.intervals.begin(), station.second.intervals.end());

    // SOLUTION/ESTIMATE
    // -----------------
    for(std::string &line : sinex.findBlock("SOLUTION/ESTIMATE")->lines)
    {
      const std::string parameterType = String::trim(line.substr(7, 6));
      const std::string name          = String::lowerCase(String::trim(line.substr(14, 4)));
      const std::string solutionId    = String::trim(line.substr(22, 4));
      const Double      value         = String::toDouble(line.substr(47, 21));

      if((stationNames.size() && std::find(stationNames.begin(), stationNames.end(), name) == stationNames.end()) ||
         (parameterType != "STAX" && parameterType != "STAY" && parameterType != "STAZ" &&
          parameterType != "VELX" && parameterType != "VELY" && parameterType != "VELZ"))
        continue;

      auto interval = std::find_if(stations[name].intervals.begin(), stations[name].intervals.end(), [&](const Interval &i) {return i.solutionId == solutionId;});
      if(interval == stations[name].intervals.end())
        throw(Exception(name+": interval for solutionId "+solutionId+" not found"));

      interval->referenceTime  = Sinex::str2time(line, 27, FALSE);
      if(     parameterType == "STAX") interval->position.x() = value;
      else if(parameterType == "STAY") interval->position.y() = value;
      else if(parameterType == "STAZ") interval->position.z() = value;
      else if(parameterType == "VELX") interval->velocity.x() = value;
      else if(parameterType == "VELY") interval->velocity.y() = value;
      else if(parameterType == "VELZ") interval->velocity.z() = value;
    }

    // discontinuities
    // ---------------
    VariableList fileNameVariableList;
    addVariable(variableLoopStation, fileNameVariableList);
    if(!fileNameDiscontinuities.empty())
    {
      logStatus<<"read discontinuity files <"<<fileNameDiscontinuities<<">"<<Log::endl;
      for(const auto &name : stationNames)
      {
        fileNameVariableList[variableLoopStation]->setValue(name);
        try
        {
          if(stations.find(name) != stations.end())
            stations[name].discontinuities = InstrumentFile::read(fileNameDiscontinuities(fileNameVariableList)).times();
        }
        catch(...)
        {
          logWarning<<name<<": discontinuity file not found"<<Log::endl;
        }
      }
    }

    if(stations.size())
    {
      logStatus<<"write instrument files <"<<fileNameInstrument<<">"<<Log::endl;
      std::vector<Time> times = timeSeries->times();
      for(auto &&station : stations)
      {
        Vector3dArc arc;
        if(!times.size() && station.second.intervals.size()) // no times given ==> use reference time of last interval
          arc = station.second.arc({station.second.intervals.back().referenceTime}, extrapolateForward, extrapolateBackward);
        else
          arc = station.second.arc(times, extrapolateForward, extrapolateBackward);
        fileNameVariableList[variableLoopStation]->setValue(station.first);
        InstrumentFile::write(fileNameInstrument(fileNameVariableList), arc);
      }
    }
    else
      logInfo<<"no stations found"<<Log::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector3dArc Sinex2StationPosition::Station::arc(const std::vector<Time> &times, Bool extrapolateForward, Bool extrapolateBackward)
{
  try
  {
    Vector3dArc arc;
    if(!intervals.size() || !times.size())
      return arc;
    if(!discontinuities.size())
      discontinuities = {Time(), date2time(2500, 1, 1)};

    // expand first/last interval in case of extrapolation so a simple interval search can be done afterwards for each epoch
    if(extrapolateBackward && times.front() < intervals.front().timeStart)
      intervals.front().timeStart = std::max(times.front(), *std::upper_bound(discontinuities.begin(), discontinuities.end(), intervals.front().timeStart));
    if(extrapolateForward && times.back() >= intervals.back().timeEnd)
      intervals.back().timeEnd = std::min(times.back()+seconds2time(1), *std::lower_bound(discontinuities.begin(), discontinuities.end(), intervals.back().timeEnd));

    UInt idInterval = 0;
    for(const auto &time : times)
    {
      while(idInterval+1 < intervals.size() && time >= intervals.at(idInterval+1).timeStart)
        idInterval++;

      if(time < intervals.at(idInterval).timeStart || time >= intervals.at(idInterval).timeEnd)
        continue; // epoch is outside interval

      Vector3dEpoch epoch;
      epoch.time = time;
      epoch.vector3d = intervals.at(idInterval).position + intervals.at(idInterval).velocity * (time-intervals.at(idInterval).referenceTime).mjd()/365.25;
      arc.push_back(epoch);
    }

    return arc;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
