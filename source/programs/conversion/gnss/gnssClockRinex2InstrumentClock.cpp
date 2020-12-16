/***********************************************/
/**
* @file gnssClockRinex2InstrumentClock.cpp
*
* @brief Convert GNSS clock RINEX files to single value instrument files for satellites or stations.
*
* @author Sebastian Strasser
* @date 2016-07-26
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring =
R"(
This program converts clocks from the \href{https://files.igs.org/pub/data/format/rinex_clock304.txt}{IGS clock RINEX format},
which contains the clocks of all satellites and stations in a single file,
into an \file{instrument file (MISCVALUE)}{instrument} for each \config{identifier}
(satellite and/or station).
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "inputOutput/file.h"
#include "files/fileInstrument.h"
#include "classes/timeSeries/timeSeries.h"

/***** CLASS ***********************************/

/** @brief Convert GNSS clock RINEX files to single value instrument files for satellites or stations.
* @ingroup programsConversionGroup */
class GnssClockRinex2InstrumentClock
{
  void readFile(const FileName &fileName, std::vector<std::string> identifier, std::vector<std::vector<Time>> &times, std::vector<std::vector<Double>> &clock) const;

public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(GnssClockRinex2InstrumentClock, SINGLEPROCESS, "Convert GNSS clock RINEX files to single value instrument files for satellites or stations.", Conversion, Gnss, Instrument)
GROOPS_RENAMED_PROGRAM(Igs2InstrumentClock, GnssClockRinex2InstrumentClock, date2time(2020, 8, 19))

/***********************************************/

void GnssClockRinex2InstrumentClock::run(Config &config)
{
  try
  {
    FileName                 fileNameOutInstrument;
    std::vector<FileName>    fileNameInClockRinex;
    TimeSeriesPtr            timesIntervalPtr;
    std::vector<std::string> identifier;
    UInt                     minEpochsPerInterval;

    renameDeprecatedConfig(config, "inputfileIgsClock", "inputfileClockRinex", date2time(2020, 8, 19));

    readConfig(config, "outputfileInstrument", fileNameOutInstrument, Config::MUSTSET,  "",  "identifier is appended to each file");
    readConfig(config, "inputfileClockRinex",  fileNameInClockRinex,  Config::MUSTSET,  "",  "");
    readConfig(config, "identifier",           identifier,            Config::MUSTSET,  "",  "satellite or station identifier, e.g. G23 or alic");
    readConfig(config, "intervals",            timesIntervalPtr,      Config::MUSTSET,  "",  "");
    readConfig(config, "minEpochsPerInterval", minEpochsPerInterval,  Config::DEFAULT,  "2", "minimum number of epochs in an interval");
    if(isCreateSchema(config)) return;

    // ======================================================

    // Read files
    // ----------
    std::vector<std::vector<Time>>   times(identifier.size());
    std::vector<std::vector<Double>> clocks(identifier.size());

    for(UInt i=0; i<fileNameInClockRinex.size(); i++)
    {
      logStatus<<"Read file <"<<fileNameInClockRinex.at(i)<<">"<<Log::endl;
      readFile(fileNameInClockRinex.at(i), identifier, times, clocks);
    }

    // Conversion
    // ----------
    const std::vector<Time> timesInterval = timesIntervalPtr->times();
    const UInt arcCount = timesInterval.size()-1;
    for(UInt idIdentifier = 0; idIdentifier < identifier.size(); idIdentifier++)
    {
      std::vector<MiscValueArc> arcs;
      UInt idx=0;

      for(UInt idArc=0; idArc<arcCount; idArc++)
      {
        // find starting epoch of arc
        while((idx<times.at(idIdentifier).size()) && (times.at(idIdentifier).at(idx)<timesInterval.at(idArc)))
          idx++;

        // fill arc with epochs
        MiscValueArc arc;
        for(UInt idEpoch=idx; (idEpoch<times.at(idIdentifier).size()) && (times.at(idIdentifier).at(idEpoch)<=timesInterval.at(idArc+1)); idEpoch++)
          if(clocks.at(idIdentifier).at(idEpoch) != 0)
          {
            MiscValueEpoch epoch;
            epoch.time = times.at(idIdentifier).at(idEpoch);
            Vector x(1);
            x(0) = clocks.at(idIdentifier).at(idEpoch);
            epoch.setData(x);
            arc.push_back(epoch);
          }

        if(arc.size()>=minEpochsPerInterval)
          arcs.push_back(arc);
        else if(arc.size())
          logWarning<<identifier.at(idIdentifier)<<": "<<timesInterval.at(idArc).dateTimeStr()<<" not enough epochs in arc ("<<arc.size()<<" < "<<minEpochsPerInterval<<")"<<Log::endl;
        else
          logWarning<<identifier.at(idIdentifier)<<": "<<timesInterval.at(idArc).dateTimeStr()<<" empty arc"<<Log::endl;
      }

      // write instrument file
      if(arcs.size())
      {
        logStatus<<"Write clocks to file <"<<fileNameOutInstrument.appendBaseName("."+identifier.at(idIdentifier))<<">"<<Log::endl;
        InstrumentFile::write(fileNameOutInstrument.appendBaseName("."+identifier.at(idIdentifier)), arcs);
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssClockRinex2InstrumentClock::readFile(const FileName &fileName, std::vector<std::string> identifier, std::vector<std::vector<Time>> &times, std::vector<std::vector<Double>> &clock) const
{
  try
  {
    InFile file(fileName);

    // read header
    Double fileVersion = 2.;
    std::string line;
    while(std::getline(file, line))
    {
      if(line.find("RINEX VERSION / TYPE") != std::string::npos)
        fileVersion = String::toDouble(line.substr(0,20));

      if(line.find("END OF HEADER") != std::string::npos)
        break;
    }

    for(auto &&id : identifier)
      std::transform(id.begin(), id.end(), id.begin(), ::toupper);

    // read data
    // ---------
    UInt v3Offset = (fileVersion >= 3.04 ? 5 : 0);
    while(std::getline(file, line))
    {
      std::string lineID = line.substr(0,7+v3Offset);
      lineID.erase(std::find_if(lineID.rbegin(), lineID.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), lineID.end()); // right trim spaces
      std::transform(lineID.begin(), lineID.end(), lineID.begin(), ::toupper);

      auto iter = std::find_if(identifier.begin(), identifier.end(), [lineID](std::string id)
      {
        return (lineID=="AS "+id || lineID=="AR "+id); }
      );
      if(iter != identifier.end())
      {
        UInt   year  = String::toInt(line.substr(v3Offset+8, 4));
        UInt   month = String::toInt(line.substr(v3Offset+13, 2));
        UInt   day   = String::toInt(line.substr(v3Offset+16, 2));
        UInt   hour  = String::toInt(line.substr(v3Offset+19, 2));
        UInt   min   = String::toInt(line.substr(v3Offset+22, 2));
        Double sec   = String::toDouble(line.substr(v3Offset+24, 10));
        Double clk   = String::toDouble(line.substr(v3Offset+40, 19));
        Time   time  = date2time(year,month,day,hour,min,sec);

        UInt idx = static_cast<UInt>(std::distance(identifier.begin(), iter));
        if(times.at(idx).size() && (time<times.at(idx).back()))
          throw(Exception("epochs not in increasing order"));

        if((times.at(idx).size()==0) || (time>times.at(idx).back()))
        {
          times.at(idx).push_back(time);
          clock.at(idx).push_back(clk);
        }
      }
    } // for(;;)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
