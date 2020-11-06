/***********************************************/
/**
* @file gnssClock2ClockRinex.cpp
*
* @brief Converts GNSS clocks from GROOPS format to IGS clock RINEX format.
*
* @author Sebastian Strasser
* @date 2019-01-22
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Converts GNSS clocks from GROOPS format to \href{ftp://igs.org/pub/data/format/rinex_clock304.txt}{IGS clock RINEX format}.
Clocks can be provided via \config{satelliteData} and/or \config{stationData}.
Observed signal types are inferred from \configFile{inputfileSignalBias}{gnssSignalBias}.
Satellites/stations used as clock references can be provided via \config{referenceClock}.

See IGS clock RINEX format description for further details on header information.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileGnssSignalBias.h"
#include "files/fileGnssStationInfo.h"
#include "files/fileInstrument.h"
#include <chrono>

/***** CLASS ***********************************/

/** @brief Converts GNSS clocks from GROOPS format to IGS clock RINEX format.
* @ingroup programsConversionGroup */
class GnssClock2ClockRinex
{
public:
  class Data
  {
  public:
    FileName inNameClock;
    std::string identifier;
    MiscValueArc arc;
    std::vector<Time> times;
    UInt idEpoch;

    Data() : idEpoch(0) {}
  };

  class SatelliteData : public Data
  {
  public:
    FileName inNameBias;

    SatelliteData() : Data() {}
  };

  class StationData : public Data
  {
  public:
    FileName inNamePosition, inNameStationInfo;
    Vector3d position;
    std::string markerNumber;

    StationData() : Data() {}
  };

  typedef std::shared_ptr<Data> DataPtr;

private:
  void readSatelliteData(std::vector<SatelliteData> &data, std::map<Char, std::set<std::string>> &system2ObsTypes) const;
  void readStationData(std::vector<StationData> &data) const;
  void writeEpoch(const Time &time, Char type, std::vector<DataPtr> &data, OutFile &file) const;

public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(GnssClock2ClockRinex, SINGLEPROCESS, "Converts GNSS clocks from GROOPS format to IGS clock RINEX format.", Conversion, Gnss, Instrument)

/***********************************************/

template<> Bool readConfig(Config &config, const std::string &name, GnssClock2ClockRinex::SatelliteData &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  if(!readConfigSequence(config, name, mustSet, defaultValue, annotation))
    return FALSE;
  readConfig(config, "inputfileClock",       var.inNameClock, Config::MUSTSET, "", "clock instrument file");
  readConfig(config, "inputfileSignalBias",  var.inNameBias,  Config::MUSTSET, "", "signal bias file");
  readConfig(config, "identifier",           var.identifier,  Config::MUSTSET, "", "PRN (e.g. G23)");
  endSequence(config);
  return TRUE;
}

/***********************************************/

template<> Bool readConfig(Config &config, const std::string &name, GnssClock2ClockRinex::StationData &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  if(!readConfigSequence(config, name, mustSet, defaultValue, annotation))
    return FALSE;
  readConfig(config, "inputfileClock",       var.inNameClock,       Config::MUSTSET, "", "clock instrument file");
  readConfig(config, "inputfilePosition",    var.inNamePosition,    Config::MUSTSET, "", "station position file");
  readConfig(config, "inputfileStationInfo", var.inNameStationInfo, Config::MUSTSET, "", "station info file");
  readConfig(config, "identifier",           var.identifier,        Config::MUSTSET, "", "station name (e.g. wtzz)");
  endSequence(config);
  return TRUE;
}

/***********************************************/

void GnssClock2ClockRinex::readSatelliteData(std::vector<SatelliteData> &data, std::map<Char, std::set<std::string>> &system2ObsTypes) const
{
  try
  {
    auto iter = data.begin();
    while(iter != data.end())
    {
      GnssSignalBias bias;
      try
      {
        iter->arc = InstrumentFile::read(iter->inNameClock);
        readFileGnssSignalBias(iter->inNameBias, bias);
      }
      catch(std::exception &e)
      {
        logWarning << e.what() << " continue..." << Log::endl;
        iter = data.erase(iter);
        continue;
      }

      iter->times = iter->arc.times();
      iter->identifier.resize(4, ' ');

      // Code biases
      for(const auto &type : bias.type)
        if(type == GnssType::RANGE)
          system2ObsTypes[type.str().at(3)].insert(type.str().substr(0,3));

      // Phase biases (for wildcard types, i.e. L1*, write one line per matching code bias type with the same phase bias; example: C1C, C1W, L1* ==> C1C, C1W, L1C, L1W)
      for(UInt i = 0; i < bias.type.size(); i++)
        if(bias.type.at(i) == GnssType::PHASE)
          for(UInt j = 0; j < bias.type.size(); j++)
            if(bias.type.at(j) == (bias.type.at(i) & ~GnssType::TYPE) + GnssType::RANGE)
              system2ObsTypes[bias.type.at(i).str().at(3)].insert(((bias.type.at(j) & ~GnssType::TYPE) + GnssType::PHASE).str().substr(0,3));
      iter++;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssClock2ClockRinex::readStationData(std::vector<StationData> &data) const
{
  try
  {
    auto iter = data.begin();
    while(iter != data.end())
    {
      GnssStationInfo info;
      try
      {
        iter->arc = InstrumentFile::read(iter->inNameClock);
        iter->position = Vector3d(InstrumentFile::read(iter->inNamePosition).at(0).data());
        readFileGnssStationInfo(iter->inNameStationInfo, info);
      }
      catch(std::exception &e)
      {
        logWarning << e.what() << " continue..." << Log::endl;
        iter = data.erase(iter);
        continue;
      }

      iter->times = iter->arc.times();
      iter->identifier.resize(4, ' ');
      std::transform(iter->identifier.begin(), iter->identifier.end(), iter->identifier.begin(), ::toupper);
      iter->markerNumber = info.markerNumber;
      iter->markerNumber.resize(20, ' ');

      iter++;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}


/***********************************************/

void GnssClock2ClockRinex::writeEpoch(const Time &time, Char type, std::vector<DataPtr> &data, OutFile &file) const
{
  try
  {
    for(auto &&d : data)
    {
      if(d->idEpoch >= d->times.size() || d->times.at(d->idEpoch) != time)
        continue; // no data for this epoch

      file << 'A' << type << ' ' << d->identifier << ' '<< time%"%y %m %d %H %M %09.6S"s << 1%"% 3i"s << "   " << d->arc.at(d->idEpoch).value%"%19.12e"s << std::endl;
      d->idEpoch++;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssClock2ClockRinex::run(Config &config)
{
  try
  {
    FileName fileNameOut;
    std::vector<SatelliteData> satelliteData;
    std::vector<StationData> stationData;
    std::vector<std::string> comment, referenceClocks;
    std::string program, institution, dcb, pcv, analysisCenter, referenceFrame;

    readConfig(config, "outputfileClockRinex",  fileNameOut,     Config::MUSTSET,   "", "");
    readConfig(config, "satelliteData",         satelliteData,   Config::OPTIONAL,  "", "one element per satellite");
    readConfig(config, "stationData",           stationData,     Config::OPTIONAL,  "", "one element per station");
    readConfig(config, "comment",               comment,         Config::OPTIONAL,  "", "comment in header");
    readConfig(config, "program",               program,         Config::MUSTSET,   "GROOPS", "name of program (for first line)");
    readConfig(config, "institution",           institution,     Config::MUSTSET,   "TUG (TU Graz)", "name of agency (for first line)");
    readConfig(config, "analysisCenter",        analysisCenter,  Config::MUSTSET,   "TUG  Graz University of Technology (TU Graz), Austria", "name of analysis center");
    readConfig(config, "differentialCodeBias",  dcb,             Config::OPTIONAL,  "GROOPS            OSB co-estimated in LSA (see bias file)", "program and source for applied differential code bias");
    readConfig(config, "phaseCenterVariations", pcv,             Config::MUSTSET,   "GROOPS            igs14.atx @ ftp.igs.org ", "program and source for applied phase center variations");
    readConfig(config, "referenceClock",        referenceClocks, Config::MUSTSET,   "", "identifier of reference satellite/station");
    readConfig(config, "referenceFrame",        referenceFrame,  Config::MUSTSET,   "IGS14", "terrestrial reference frame for the stations");
    if(isCreateSchema(config)) return;

    logStatus << "read data files" << Log::endl;
    std::map<Char, std::set<std::string>> system2ObsTypes;
    readSatelliteData(satelliteData, system2ObsTypes);
    readStationData(stationData);

    if(!satelliteData.size() && !stationData.size())
      throw(Exception("no satellite or station data provided"));

    logStatus << "write clock RINEX file <" << fileNameOut << ">" << Log::endl;
    OutFile file(fileNameOut);

    // Header: RINEX VERSION / TYPE
    file << 3.02%"%9.2f"s << std::string(11, ' ') << 'C' << std::string(19, ' ') << (system2ObsTypes.size() == 1 ? system2ObsTypes.begin()->first : 'M')
         << std::string(19, ' ') << "RINEX VERSION / TYPE" << std::endl;

    // Header: PGM / RUN BY / DATE
    auto systemTime = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    program.resize(20, ' ');
    institution.resize(20, ' ');
    file << program << institution << std::put_time(gmtime(&systemTime), "%Y%m%d %H%M%S") << " UTC PGM / RUN BY / DATE" << std::endl;

    // Header: COMMENT
    for(auto &c : comment)
    {
      c.resize(60, ' ');
      file << c << "COMMENT" << std::endl;
    }

    // Header: SYS / # / OBS TYPES
    for(const auto &system : system2ObsTypes)
    {
      file << system.first << "  " << system.second.size()%"% 3i"s;
      std::vector<std::string> types(system.second.begin(), system.second.end());
      types.resize(static_cast<UInt>(std::ceil(types.size()/13.)*13), "   ");
      for(UInt i = 0; i < types.size(); i++)
      {
        file << " " << types.at(i);
        if(i > 0 && (i+1)%13 == 0)
          file << "  SYS / # / OBS TYPES" << std::endl << (i+1 < types.size() ? "      " : "");
      }
    }

    // Header: TIME SYSTEM ID
    file << std::string(3, ' ') << "GPS" << std::string(54, ' ') << "TIME SYSTEM ID" << std::endl;

    // Header: SYS / DCBS APPLIED
    dcb.resize(58, ' ');
    for(const auto &system : system2ObsTypes)
      file << system.first << " " << dcb << "SYS / DCBS APPLIED" << std::endl;

    // Header: SYS / PCVS APPLIED
    pcv.resize(58, ' ');
    for(const auto &system : system2ObsTypes)
      file << system.first << " " << pcv << "SYS / PCVS APPLIED" << std::endl;

    // Header: # / TYPES OF DATA
    std::vector<std::string> dataTypes;
    if(satelliteData.size())
      dataTypes.push_back("AS");
    if(stationData.size())
      dataTypes.push_back("AR");
    file << dataTypes.size()%"% 6i"s;
    dataTypes.resize(5, "  ");
    for(const auto &type : dataTypes)
      file << std::string(4, ' ') << type;
    file << std::string(24, ' ') << "# / TYPES OF DATA" << std::endl;

    // Header: ANALYSIS CENTER
    analysisCenter.resize(60, ' ');
    file << analysisCenter << "ANALYSIS CENTER" << std::endl;

    // Header: # OF CLK REF
    file << referenceClocks.size()%"% 6i"s << std::string(54, ' ') << "# OF CLK REF " << std::endl;

    // Header: ANALYSIS CLK REF
    for(auto ref : referenceClocks)
    {
      ref.resize(60, ' ');
      file << ref << "ANALYSIS CLK REF" << std::endl;
    }

    // Header: # OF SOLN STA / TRF
    if(stationData.size())
    {
      referenceFrame.resize(50, ' ');
      file << stationData.size()%"% 6i"s << std::string(4, ' ') << referenceFrame << "# OF SOLN STA / TRF" << std::endl;
    }

    // Header: SOLN STA NAME / NUM
    if(stationData.size())
      for(const auto &station : stationData)
        file << station.identifier << " " << station.markerNumber << (station.position.x()*1e3)%"%11.0f"s << " "
             << (station.position.y()*1e3)%"%11.0f"s << " " << (station.position.z()*1e3)%"%11.0f"s << "SOLN STA NAME / NUM" << std::endl;

    // Header: # OF SOLN SATS
    if(satelliteData.size())
      file << satelliteData.size()%"% 6i"s << std::string(54, ' ') << "# OF SOLN SATS " << std::endl;

    // Header: PRN LIST
    const UInt count = std::ceil(satelliteData.size()/15.)*15;
    for(UInt i = 0; i < count; i++)
    {
      file << (i < satelliteData.size() ? satelliteData.at(i).identifier : std::string(4, ' '));
      if(i > 0 && (i+1)%15 == 0)
        file << "PRN LIST" << std::endl;
    }

    file << std::string(60, ' ') << "END OF HEADER" << std::endl;

    // Clock data
    std::set<Time> times; // unique set of times covering all satellites and stations
    for(auto &&data : satelliteData)
      std::copy(data.times.begin(), data.times.end(), std::inserter(times, times.end()));
    for(auto &&data : stationData)
      std::copy(data.times.begin(), data.times.end(), std::inserter(times, times.end()));

    UInt i = 0;
    std::vector<DataPtr> satelliteDataPtrs, stationDataPtrs;
    for(auto & data : satelliteData)
      satelliteDataPtrs.push_back(std::make_shared<Data>(data));
    for(auto & data : stationData)
      stationDataPtrs.push_back(std::make_shared<Data>(data));
    logTimerStart
    for(const auto &time : times)
    {
      logTimerLoop(++i, times.size())
      writeEpoch(time, 'S', satelliteDataPtrs, file);
      writeEpoch(time, 'R', stationDataPtrs,   file);
    }
    logTimerLoopEnd(times.size());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
