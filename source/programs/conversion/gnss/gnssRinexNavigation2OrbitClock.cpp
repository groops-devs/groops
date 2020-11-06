/***********************************************/
/**
* @file gnssRinexNavigation2OrbitClock.cpp
*
* @brief Convert RINEX navigation file (e.g. broadcast ephemeris) to orbit and clock files.
*
* @author Sebastian Strasser
* @date 2017-02-28
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Evaluates orbit and clock parameters from \href{ftp://ftp.igs.org/pub/data/format/rinex304.pdf}{RINEX navigation file} \config{inputfileRinex}
at epochs given by \configClass{timeSeries}{timeSeriesType} and writes them to \configFile{outputfileOrbit}{instrument} and
\configFile{outputfileClock}{instrument}, respectively.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "inputOutput/file.h"
#include "files/fileInstrument.h"
#include "classes/timeSeries/timeSeries.h"

/***** CLASS ***********************************/

/** @brief Convert RINEX navigation file (e.g. broadcast ephemeris) to orbit and clock files.
* @ingroup programsConversionGroup */
class GnssRinexNavigation2OrbitClock
{
  class Satellite
  {
  public:
    std::vector<Time>   epoch;
    std::vector<Vector> clockParam;
    std::vector<Matrix> orbitParam;
  };

  Double rinexVersion;
  Char   system;
  Vector ionAlpha, ionBeta, deltaUTC;
  Double leapSeconds;
  std::map<GnssType, Satellite> satellites;
  std::vector<Time> times;

  void readHeader(InFile &file, UInt lineCount=MAX_UINT);
  void readData(InFile &file);
  Bool getLine(InFile &file, std::string &line, std::string &label) const;
  Bool testLabel(const std::string &labelInLine, const std::string &label, Bool optional=TRUE) const;

public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(GnssRinexNavigation2OrbitClock, SINGLEPROCESS, "Convert RINEX navigation file (e.g. broadcast ephemeris) to orbit and clock files.", Conversion, Gnss, Instrument)

/***********************************************/

void GnssRinexNavigation2OrbitClock::run(Config &config)
{
  try
  {
    FileName outNameOrbit, outNameClock, inNameRinex;
    TimeSeriesPtr timeSeriesPtr;
    std::vector<std::string> usePrn;

    readConfig(config, "outputfileOrbit", outNameOrbit,  Config::OPTIONAL, "", "PRN is appended to file name");
    readConfig(config, "outputfileClock", outNameClock,  Config::OPTIONAL, "", "PRN is appended to file name");
    readConfig(config, "inputfileRinex",  inNameRinex,   Config::MUSTSET,  "", "RINEX navigation file");
    readConfig(config, "timeSeries",      timeSeriesPtr, Config::MUSTSET,  "", "orbit and clock evaluation epochs");
    readConfig(config, "usePrn",          usePrn,        Config::OPTIONAL, "", "only export these PRNs instead of all");
    if(isCreateSchema(config)) return;

    times = timeSeriesPtr->times();

    // read file
    logStatus << "read RINEX file <" << inNameRinex << ">" << Log::endl;
    InFile file(inNameRinex);
    readHeader(file);
    readData(file);
    file.close();

    if(!outNameOrbit.empty())
    {
      throw(Exception("orbit evaluation not implemented yet"));
    }

    if(!outNameClock.empty())
    {
      logStatus << "write clock files <" << outNameClock.appendBaseName(".***") << ">" << Log::endl;
      for(const auto &satellite : satellites)
      {
        if(!satellite.second.epoch.size() || (usePrn.size() && std::find(usePrn.begin(), usePrn.end(), satellite.first.str().substr(3,3)) == usePrn.end()))
          continue;

        // fill arc with epochs
        UInt idx = 0;
        MiscValueArc arc;
        for(UInt idEpoch = 0; idEpoch < times.size(); idEpoch++)
        {
          while(idx+1 < satellite.second.epoch.size() && times.at(idEpoch) >= satellite.second.epoch.at(idx+1))
            idx++;

          MiscValueEpoch epoch;
          epoch.time = times.at(idEpoch);
          const Double dt = (times.at(idEpoch) - satellite.second.epoch.at(idx)).seconds();
          for(UInt i = 0; i < satellite.second.clockParam.at(idx).size(); i++)
            epoch.value += satellite.second.clockParam.at(idx)(i)*std::pow(dt, i);  // a_i*dt^i
          arc.push_back(epoch);
        }

        // write clock file
        std::list<Arc> arcList;
        arcList.push_back(arc);
        InstrumentFile::write(outNameClock.appendBaseName('.'+satellite.first.str().substr(3,3)), arcList);
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssRinexNavigation2OrbitClock::readHeader(InFile &file, UInt lineCount)
{
  try
  {
    std::string line, label;
    getLine(file, line, label);
    testLabel(label, "RINEX VERSION / TYPE", FALSE);
    rinexVersion    = String::toDouble(line.substr(0, 9));
    if(rinexVersion<2)
      logWarning << "old RINEX version: " << rinexVersion << Log::endl;
    if(line.at(20) != 'N' && !(rinexVersion < 3 && line.at(20) == 'G'))
      throw(Exception("File must contain Navigation Data"));
    system = rinexVersion < 3 ? (line.at(20) == 'N' ? 'G' : 'R') : line.at(40);

    for(UInt idLine=0; idLine<lineCount; idLine++)
    {
      if(!getLine(file, line, label))
        throw(Exception("error while reading RINEX header"));
      if(std::all_of(line.begin(), line.end(), isspace))
      {
        if(rinexVersion < 2)
          break;
        else
          continue;
      }
      // ====================================
      if(testLabel(label, "END OF HEADER"))
        break;
      // ====================================
      else if(testLabel(label, "PGM / RUN BY / DATE"))
      {
      }
      // ====================================
      else if(testLabel(label, "COMMENT"))
      {
      }
      // ====================================
      else if(testLabel(label, "ION ALPHA"))
      {
        ionAlpha = Vector(4);
        for(UInt i = 0; i < ionAlpha.size(); i++)
          ionAlpha(i) = String::toDouble(line.substr(2+i*12, 12));
      }
      // ====================================
      else if(testLabel(label, "ION BETA"))
      {
        ionBeta = Vector(4);
        for(UInt i = 0; i < ionBeta.size(); i++)
          ionBeta(i) = String::toDouble(line.substr(2+i*12, 12));
      }
      // ====================================
      else if(testLabel(label, "DELTA-UTC: A0,A1,T,W"))
      {
        deltaUTC = Vector(4);
        deltaUTC(0) = String::toDouble(line.substr(3, 19));
        deltaUTC(1) = String::toDouble(line.substr(22, 19));
        deltaUTC(2) = String::toInt(line.substr(41, 9));
        deltaUTC(3) = String::toInt(line.substr(50, 9));
      }
      // ====================================
      else if(testLabel(label, "LEAP SECONDS"))
      {
        leapSeconds = String::toInt(line.substr(0, 6));
      }
      // ====================================
      else if(testLabel(label, "CORR TO SYSTEM TIME"))
      {
      }
      // ====================================
      else if(testLabel(label, "TIME SYSTEM CORR"))
      {
      }
      // ====================================
      else if(testLabel(label, "IONOSPHERIC CORR"))
      {
      }
      // ====================================
      else
      {
        logWarning<<"Unknown header label:"<<Log::endl;
        logWarning<<"'"<<line<<"'"<<Log::endl;
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssRinexNavigation2OrbitClock::readData(InFile &file)
{
  std::string line, label;

  try
  {
    while(getLine(file, line, label))
    {
      std::string prnStr = rinexVersion < 3 ? system+line.substr(0,2): line.substr(0,3);
      if(prnStr.at(1) == ' ') prnStr.at(1) = '0';
      GnssType prn("***" + prnStr);

      const UInt lineCount = (prn == GnssType::GLONASS || prn == GnssType::SBAS) ? 3 : (rinexVersion < 2 ? 6 : 7);

      // epoch
      Int year   = String::toInt(line.substr(rinexVersion < 3 ?  3 :  4, rinexVersion < 3 ? 2 : 4));
      Int month  = String::toInt(line.substr(rinexVersion < 3 ?  6 :  9, 2));
      Int day    = String::toInt(line.substr(rinexVersion < 3 ?  9 : 12, 2));
      Int hour   = String::toInt(line.substr(rinexVersion < 3 ? 12 : 15, 2));
      Int minute = String::toInt(line.substr(rinexVersion < 3 ? 15 : 18, 2));
      Double sec = String::toDouble(line.substr(rinexVersion < 3 ? 17 : 21, rinexVersion < 3 ? 5 : 2));
      if(rinexVersion < 3)
        year += ((year<=80) ? 2000 : 1900);
      const Time time = date2time(year, month, day, hour, minute, sec);

      if(satellites[prn].epoch.size() && time <= satellites[prn].epoch.back())
      {
        //logWarning << prn.str().substr(3,3) << ": duplicate entry at " << time.dateTimeStr() << Log::endl;
        for(UInt i = 0; i < lineCount; i++)
          getLine(file, line, label);
        continue;
      }
      satellites[prn].epoch.push_back(time);

      // clock polynomial
      Vector clockParam((prn != GnssType::GLONASS && prn != GnssType::SBAS) ? 3 : 2);
      for(UInt i = 0; i < clockParam.size(); i++)
        clockParam(i) = String::toDouble(line.substr((rinexVersion < 3 ? 22 : 23)+i*19, 19));
      satellites[prn].clockParam.push_back(clockParam);

      // orbit parameters
      Matrix orbitParam(lineCount,4);
      for(UInt i = 0; i < lineCount; i++)
      {
        getLine(file, line, label);
        for(UInt j = 0; j < 4; j++)
          orbitParam(i,j) = String::toDouble(line.substr((rinexVersion < 3 ? 3 : 4)+j*19, 19));
      }
      satellites[prn].orbitParam.push_back(orbitParam);
    }
  }
  catch(std::exception &e)
  {
    logError<<"'"<<line<<"'"<<Log::endl;
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool GnssRinexNavigation2OrbitClock::getLine(InFile &file, std::string &line, std::string &label) const
{
  try
  {
    getline(file, line);
    if(line.back() == '\r')
      line.pop_back();
    if(line.size()<80)
      line.resize(80,' ');
    label = line.substr(60,20);
    return file.good();
  }
  catch(...)
  {
    line.clear();
    line.resize(80,' ');
    label = line.substr(60,20);
    return FALSE;
  }
}

/***********************************************/

Bool GnssRinexNavigation2OrbitClock::testLabel(const std::string &labelInLine, const std::string &label, Bool optional) const
{
  if(labelInLine.find(label)!=std::string::npos)
    return TRUE;
  if(optional)
    return FALSE;
  throw(Exception(std::string("In Line '")+labelInLine+"' label '"+label+"' expected\n"));
}

/***********************************************/
