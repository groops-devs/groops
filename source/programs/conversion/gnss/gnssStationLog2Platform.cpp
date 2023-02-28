/***********************************************/
/**
* @file gnssStationLog2Platform.cpp
*
* @brief GNSS analysis.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2012-04-28
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Converts \href{https://files.igs.org/pub/station/general/blank.log}{IGS station log format} to \configFile{outputfileStationPlatform}{platform}.

If \configFile{inputfileAntennaDefinition}{gnssAntennaDefinition} is provided, station log data is cross-checked with the given antenna definitions.
Cross-checking station log data with a \href{https://www.iers.org/IERS/EN/Organization/AnalysisCoordinator/SinexFormat/sinex.html}{SINEX file} is also
possible by providing \config{inputfileSinex}. Any failed checks result in warnings in the output log.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "inputOutput/file.h"
#include "inputOutput/fileSinex.h"
#include "files/filePlatform.h"

/***** CLASS ***********************************/

/** @brief GNSS analysis.
* @ingroup programsConversionGroup */
class GnssStationLog2Platform
{
  Platform readFile(const FileName &fileName);
  void     newAntenna (Platform &platform, PlatformGnssAntenna &antenna) const;
  void     newReceiver(Platform &platform, PlatformGnssReceiver &receiver) const;
  Bool     readDouble(const std::string &line, Double &x) const;
  Bool     readString(const std::string &line, std::string &x) const;
  Bool     readTime  (const std::string &line, Time &x) const;
  void     checkSinexFile(const FileName &fileName, Platform &platform) const;

public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GnssStationLog2Platform, SINGLEPROCESS, "GNSS analysis", Conversion, Gnss)
// GROOPS_RENAMED_PROGRAM(GnssStationLog2StationInfo, GnssStationLog2Platform, date2time(2023, 1, 4))

/***********************************************/

void GnssStationLog2Platform::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameStationLog, fileNameAntenna, fileNameStationPlatform, fileNameSinex;

    renameDeprecatedConfig(config, "outputfileStationInfo", "outputfileStationPlatform", date2time(2023, 1, 4));

    readConfig(config, "outputfileStationPlatform",  fileNameStationPlatform, Config::MUSTSET,  "", "");
    readConfig(config, "inputfileStationLog",        fileNameStationLog,      Config::MUSTSET,  "", "");
    readConfig(config, "inputfileAntennaDefinition", fileNameAntenna,         Config::OPTIONAL, "{groopsDataDir}/gnss/receiverStation/antennaDefinition/igs/igs20/antennaDefinition_igs20.dat", "used to check antennas");
    readConfig(config, "inputfileSinex",             fileNameSinex,           Config::OPTIONAL, "{groopsDataDir}/gnss/receiverStation/igs_with_former.snx", "used to cross-check station log with SINEX file");
    if(isCreateSchema(config)) return;

    // ==========================

    logStatus<<"read station log <"<<fileNameStationLog<<">"<<Log::endl;
    Platform platform = readFile(fileNameStationLog);

    // ==========================

    // some tests
    // ----------
    if(platform.markerName.size()!=4)
      throw(Exception(platform.markerName+"."+platform.markerNumber+": marker name should have 4 letters"));
    if(!std::any_of(platform.equipments.begin(), platform.equipments.end(), [](const auto &p) {return std::dynamic_pointer_cast<PlatformGnssReceiver>(p);}))
      throw(Exception(platform.markerName+"."+platform.markerNumber+": No receiver given"));
    if(!std::any_of(platform.equipments.begin(), platform.equipments.end(), [](const auto &p) {return std::dynamic_pointer_cast<PlatformGnssAntenna>(p);}))
      throw(Exception(platform.markerName+"."+platform.markerNumber+": No antenna given"));

    // ==========================

    // test antennas
    // -------------
    if(!fileNameAntenna.empty())
    {
      logStatus<<"read antenna definition file <"<<fileNameAntenna<<">"<<Log::endl;
      std::vector<GnssAntennaDefinitionPtr> antennaList;
      readFileGnssAntennaDefinition(fileNameAntenna, antennaList);

      platform.fillGnssAntennaDefinition(antennaList);
      // test antennas
      for(const auto &equipment : platform.equipments)
      {
        auto antenna = std::dynamic_pointer_cast<PlatformGnssAntenna>(equipment);
        if(antenna && !antenna->antennaDef)
          logWarning<<platform.markerName<<"."<<platform.markerNumber<<": no antenna definition found for "
                    <<antenna->str()<<Log::endl;
      }
    }

    // ==========================

    // check platform with sinex file
    // -----------------------------
    if(!fileNameSinex.empty())
      checkSinexFile(fileNameSinex, platform);

    // test if there is still no approx position (neither in station log nor SINEX file)
    if(platform.approxPosition.r()<6300e3)
    {
      logWarning<<platform.markerName<<"."<<platform.markerNumber<<": No approx. position given"<<Log::endl;
      throw(Exception("No approx. position given"));
    }

    // ==========================

    logStatus<<"write station info <"<<fileNameStationPlatform<<">"<<Log::endl;
    writeFilePlatform(fileNameStationPlatform, platform);

    // test platform
    try
    {
      readFilePlatform(fileNameStationPlatform, platform);
    }
    catch(std::exception &/*e*/)
    {
      logWarning<<"not readable: '"<<fileNameStationPlatform<<"'"<<Log::endl;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

Platform GnssStationLog2Platform::readFile(const FileName &fileName)
{
  std::string line;

  try
  {
    Platform  platform;
    PlatformGnssReceiver receiver;
    PlatformGnssAntenna  antenna;

    InFile file(fileName);

    enum Section {OTHER    = 999,
                  SITE     = 1,
                  LOCATION = 2,
                  RECEIVER = 3,
                  ANTENNA  = 4};
    Section section = OTHER;

    while(std::getline(file, line))
    {
      if(line.empty())
        continue;

      if(line[line.length()-1] == '\r') {
        line.erase(line.length()-1);
      }

      UInt pos = line.find_first_not_of(" \t\n\r");
      if(pos == std::string::npos || pos > 30)
        continue;
      Char firstLetter = line.at(pos);

      // section change?
      // ---------------
      if(line.at(0) == '1') // GNSS Receiver Information
      {
        section = SITE;
      }
      else if(line.at(0) == '2') // GNSS Receiver Information
      {
        section = LOCATION;
      }
      else if(line.at(0) == '3') // GNSS Receiver Information
      {
        section = RECEIVER;

        newReceiver(platform, receiver);
        receiver = PlatformGnssReceiver();
      }
      else if(line.at(0) == '4')  // GNSS Antenna Information
      {
        section = ANTENNA;

        newAntenna(platform, antenna);
        antenna = PlatformGnssAntenna();
      }
      else if(line.at(0) == '5')  // GNSS Antenna Information
      {
        break;
      }
      else if((line.at(0) != ' ')&&(line.at(0) != '\t')) // other section
      {
        section = OTHER;
      }

      // =================================================

      // skip other sections
      if(section == OTHER)
        continue;

      // =================================================

      if(section == SITE)
      {
        if(line.find("Site Identification of the GNSS Monument") != std::string::npos)
        {
        }
        else if(line.find("Four Character ID") != std::string::npos)
        {
          if(platform.markerName.empty())
          {
            readString(line, platform.markerName);
            std::transform(platform.markerName.begin(), platform.markerName.end(), platform.markerName.begin(), ::toupper);
            if(platform.markerName.length() > 4)
            {
              logWarning<<platform.markerName<<": marker name has more than 4 letters, shortened to: "<<platform.markerName.substr(0,4)<<Log::endl;
              platform.markerName = platform.markerName.substr(0,4);
            }
          }
        }
        else if(line.find("IERS DOMES Number") != std::string::npos)
        {
          if(platform.markerNumber.empty())
            readString(line, platform.markerNumber);
        }
/*        else
          logInfo<<"unknown line '"<<line<<"'"<<Log::endl;*/
      } // if(section == SITE)

      // =================================================

      if(section == LOCATION)
      {
        if(line.find("Site Location Information") != std::string::npos)
        {
        }
        else if(line.find("City or Town") != std::string::npos)
        {
          readString(line, platform.comment);
        }
        else if(line.find("Country") != std::string::npos)
        {
          std::string country;
          readString(line, country);
          platform.comment += ", " + country;
        }
        else if(line.find("X coordinate") != std::string::npos)
        {
          readDouble(line, platform.approxPosition.x());
        }
        else if(line.find("Y coordinate") != std::string::npos)
        {
          readDouble(line, platform.approxPosition.y());
        }
        else if(line.find("Z coordinate") != std::string::npos)
        {
          readDouble(line, platform.approxPosition.z());
        }
        else if(platform.approxPosition.r() < 6300e3 && line.find("Latitude") != std::string::npos)
        {
          readDouble(line, platform.approxPosition.x());
        }
        else if(platform.approxPosition.r() < 6300e3 && line.find("Longitude") != std::string::npos)
        {
          readDouble(line, platform.approxPosition.y());
        }
        else if(platform.approxPosition.r() < 6300e3 && line.find("Elevation") != std::string::npos)
        {
          readDouble(line, platform.approxPosition.z());
        }
/*        else
          logInfo<<"unknown line '"<<line<<"'"<<Log::endl;*/
      } // if(section == LOCATION)

      // =================================================

      if(section == RECEIVER)
      {
        if((line.find("Additional Information") != std::string::npos) || (firstLetter==':'))
        {
//           std::string comment;
//           readString(line, comment);
//           receiver.comment += comment;
        }
        else if(line.find("GNSS Receiver Information") != std::string::npos)
        {
        }
        else if(line.find("Receiver Type") != std::string::npos)
        {
          readString(line, receiver.name);
        }
        else if(line.find("Satellite System") != std::string::npos)
        {
        }
        else if(line.find("Serial Number") != std::string::npos)
        {
          readString(line, receiver.serial);
        }
        else if(line.find("Firmware Version") != std::string::npos)
        {
          readString(line, receiver.version);
        }
        else if(line.find("Elevation Cutoff") != std::string::npos)
        {
        }
        else if(line.find("Date Installed") != std::string::npos)
        {
          readTime(line, receiver.timeStart);
        }
        else if(line.find("Date Removed") != std::string::npos)
        {
          readTime(line, receiver.timeEnd);
        }
        else if(line.find("Temperature Stabiliz") != std::string::npos)
        {
        }
//         else
//           logWarning<<platform.markerName<<"."<<platform.markerNumber<<" unknown line in GNSS Receiver Information section '"<<line<<"'"<<Log::endl;
      } // if(section == RECEIVER)

      // =================================================

      if(section == ANTENNA)
      {
        if((line.find("Additional Information") != std::string::npos) || (firstLetter==':'))
        {
//           std::string comment;
//           readString(line, comment);
//           antenna.comment += comment;
        }
        else if(line.find("GNSS Antenna Information") != std::string::npos)
        {
        }
        else if(line.find("Antenna Type") != std::string::npos)
        {
          std::string text;
          readString(line, text);
          if(text.size() == 20) // with radome?
          {
            antenna.radome = text.substr(16, 4);
            if(antenna.radome == "NONE")
              antenna.radome.clear();
          }
          if(text.substr(0, 8) == "TRIMBLE ")
            text = text.substr(8);
          if(text == "DORNE MARGOLIN T")
            text = "AOAD/M_T";
          text = String::trim(text.substr(0, 15));
          antenna.name = text;
        }
        else if(line.find("Radome Serial Number") != std::string::npos)
        {
        }
        else if(line.find("Serial Number") != std::string::npos)
        {
          readString(line, antenna.serial);
          if((antenna.serial=="n/a")||(antenna.serial=="N/A"))
            antenna.serial.clear();
        }
        else if(line.find("Antenna Reference Point") != std::string::npos)
        {
        }
        else if(line.find("Marker->ARP Up") != std::string::npos)
        {
          readDouble(line, antenna.position.z());
        }
        else if(line.find("Marker->ARP East") != std::string::npos)
        {
          readDouble(line, antenna.position.y());
        }
        else if(line.find("Marker->ARP North") != std::string::npos)
        {
          readDouble(line, antenna.position.x());
        }
        else if(line.find("Alignment from True N") != std::string::npos || line.find("Degree Offset from North") != std::string::npos /*old format*/)
        {
          Double deg;
          if(readDouble(line, deg) && (deg!=0))
            antenna.local2antennaFrame = rotaryZ(Angle(deg*DEG2RAD));
        }
        else if(line.find("Antenna Radome Type") != std::string::npos)
        {
          std::string text;
          readString(line, text);
          if((!antenna.radome.empty()) && (text != antenna.radome))
            logWarning<<"Radome different: "<<antenna.radome<<" != "<<text<<Log::endl;
          antenna.radome = text;
          std::transform(antenna.radome.begin(), antenna.radome.end(), antenna.radome.begin(), ::toupper);
          if(antenna.radome == "NONE")
            antenna.radome.clear();
        }
        else if(line.find("Antenna Cable Type") != std::string::npos)
        {
        }
        else if(line.find("Antenna Cable Length") != std::string::npos)
        {
        }
        else if(line.find("Date Installed") != std::string::npos)
        {
          readTime(line, antenna.timeStart);
        }
        else if(line.find("Date Removed") != std::string::npos)
        {
          readTime(line, antenna.timeEnd);
        }
        else if(line.find("Antenna Height") != std::string::npos) // old format
        {
          readDouble(line, antenna.position.z());
        }
//         else
//           logWarning<<platform.markerName<<"."<<platform.markerNumber<<" unknown line in GNSS Antenna Information section '"<<line<<"'"<<Log::endl;
      } // if(section == ANTENNA)

    } // for(;;)

    // if no Cartesian coordinates are given, try to use ellipsoidal coordinates
    if(platform.approxPosition.r() > 0 && platform.approxPosition.r() < 6300e3)
    {
      Ellipsoid ellipsoid;
      platform.approxPosition = ellipsoid(Angle(platform.approxPosition.y()*DEG2RAD), Angle(platform.approxPosition.x()*DEG2RAD), platform.approxPosition.z());
    }

    // insert last receiver
    // --------------------
    newReceiver(platform, receiver);

    // insert last antenna
    // -------------------
    newAntenna(platform, antenna);

    return platform;
  }
  catch(std::exception &e)
  {
    logError<<"error: '"<<line<<"'"<<Log::endl;
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssStationLog2Platform::newReceiver(Platform &platform, PlatformGnssReceiver &receiver) const
{
  try
  {
    if(receiver.name.empty() || receiver.name == "no receiver")
      return;

    if(receiver.timeStart == Time())
      logWarning<<platform.markerName<<"."<<platform.markerNumber<<" receiver timeStart not given"<<Log::endl;

    if(receiver.timeEnd == Time())
      receiver.timeEnd = date2time(2500,1,1);
    if(receiver.timeStart >= receiver.timeEnd)
      logWarning<<platform.markerName<<"."<<platform.markerNumber<<": receiver time error (new receiver start: "
                <<receiver.timeStart.dateTimeStr()<<", new receiver end: "<<receiver.timeEnd.dateTimeStr()<<")"<<Log::endl;

    // find last receiver
    auto iter = std::find_if(platform.equipments.rbegin(), platform.equipments.rend(),
                             [](const auto x) {return x->getType() == PlatformEquipment::GNSSRECEIVER;});
    if(iter == platform.equipments.rend())
    {
      platform.equipments.push_back(std::make_shared<PlatformGnssReceiver>(receiver));
      return;
    }
    std::shared_ptr<PlatformGnssReceiver> lastReceiver = std::dynamic_pointer_cast<PlatformGnssReceiver>(*iter);

    if((lastReceiver->name    == receiver.name)   &&
       (lastReceiver->serial  == receiver.serial) &&
       (lastReceiver->version == receiver.version))
    {
      lastReceiver->timeEnd = receiver.timeEnd;
      return;
    }

    // test times
    if(lastReceiver->timeEnd == date2time(2500,1,1))
      lastReceiver->timeEnd = receiver.timeStart;

    if(lastReceiver->timeEnd > receiver.timeStart)
      logWarning<<platform.markerName<<"."<<platform.markerNumber<<" receiver time error (previous receiver end: "
               <<lastReceiver->timeEnd.dateTimeStr()<<", new receiver start: "<<receiver.timeStart.dateTimeStr()<<")"<<Log::endl;

    platform.equipments.push_back(std::make_shared<PlatformGnssReceiver>(receiver));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssStationLog2Platform::newAntenna(Platform &platform, PlatformGnssAntenna &antenna) const
{
  try
  {
    if(antenna.name.empty() || antenna.name == "no antenna")
      return;

    if(antenna.timeStart == Time())
      logWarning<<platform.markerName<<"."<<platform.markerNumber<<" antenna timeStart not given"<<Log::endl;

    if(antenna.timeEnd == Time())
      antenna.timeEnd = date2time(2500,1,1);
    if(antenna.timeStart >= antenna.timeEnd)
      logWarning<<platform.markerName<<"."<<platform.markerNumber<<" antenna time error"<<Log::endl;

    // find last antenna
    auto iter = std::find_if(platform.equipments.rbegin(), platform.equipments.rend(),
                             [](const auto x) {return x->getType() == PlatformEquipment::GNSSANTENNA;});
    if(iter == platform.equipments.rend())
    {
      platform.equipments.push_back(std::make_shared<PlatformGnssAntenna>(antenna));
      return;
    }
    auto lastAntenna = std::dynamic_pointer_cast<PlatformGnssAntenna>(*iter);

    if((lastAntenna->name    == antenna.name)   &&
       (lastAntenna->serial  == antenna.serial) &&
       (lastAntenna->radome  == antenna.radome) &&
       ((lastAntenna->position-antenna.position).r()<0.0001) &&
       quadsum(lastAntenna->local2antennaFrame.matrix() - antenna.local2antennaFrame.matrix())<0.0001)
    {
      lastAntenna->timeEnd = antenna.timeEnd;
      return;
    }

    // test times
    if(lastAntenna->timeEnd == date2time(2500,1,1))
      lastAntenna->timeEnd = antenna.timeStart - seconds2time(1.);

    if(lastAntenna->timeEnd > antenna.timeStart)
      logWarning<<platform.markerName<<"."<<platform.markerNumber<<" antenna time error"<<Log::endl;

    platform.equipments.push_back(std::make_shared<PlatformGnssAntenna>(antenna));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool GnssStationLog2Platform::readDouble(const std::string &line, Double &x) const
{
  try
  {
    UInt pos = line.find(":");
    if(pos == std::string::npos)
      throw(Exception("':' expected"));
    std::string text = line.substr(pos+1);
    std::stringstream ss(text);
    Double f = 0;
    ss>>f;
    if((!ss.good())&&(!ss.eof()))
      return FALSE;
    x = f;
    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool GnssStationLog2Platform::readString(const std::string &line, std::string &x) const
{
  try
  {
    UInt pos = line.find(":");
    if(pos == std::string::npos)
      throw(Exception("':' expected"));
    std::string text = line.substr(pos+1);
    if(text.find("(") != std::string::npos && text.find("(") == 1)
      return FALSE;
    UInt start = text.find_first_not_of(" \t\n\r");
    if(start == std::string::npos)
      return FALSE;
    UInt end = text.find_last_not_of(" \t\n\r");
    x = text.substr(start, end-start+1);
    std::replace_if(x.begin(), x.end(), [](char c){ return c < 0; }, '?'); // replace non ASCII characters
    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool GnssStationLog2Platform::readTime(const std::string &line, Time &x) const
{
  try
  {
    UInt pos = line.find(":");
    if(pos == std::string::npos)
      throw(Exception("':' expected"));
    std::string text = line.substr(pos+1);
    if(text.find("(") != std::string::npos)
      return FALSE;
    UInt start = text.find_first_not_of(" \t\n\r");
    if(start == std::string::npos)
      return FALSE;
    std::stringstream ss(text);

    UInt year, month, day;
    UInt hour = 0, min = 0;
    char c;
    ss>>year;

    if(year < 1950) // old date style: 02-OCT-2000 00:00 UT
    {
      day = year;
      if((!ss.good())&&(!ss.eof()))
        return FALSE;
      ss>>c; if(c!='-') return FALSE;
      char tmp[3];
      ss.read(tmp,3);
      std::string monthName(tmp,3);
      std::transform(monthName.begin(), monthName.end(), monthName.begin(), ::toupper);
      const std::vector<std::string> months = {"JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"};
      auto pos = std::find(months.begin(), months.end(), monthName);
      if(pos == months.end())
        return FALSE;
      month = std::distance(months.begin(), pos)+1;
      ss>>c; if(c!='-') return FALSE;
      ss>>year;
      if((!ss.good())&&(!ss.eof()))
        return FALSE;
      ss>>std::noskipws;
      ss>>c;
      if(c==' ')
      {
      ss>>hour;
      if((!ss.good())&&(!ss.eof()))
        return FALSE;
      ss>>c; if(c!=':') return FALSE; //throw(Exception(": expected"));
      ss>>min;
      if((!ss.good())&&(!ss.eof()))
        return FALSE;
      }
    }
    else // new date style: 2015-01-28T00:00Z
    {
      if((!ss.good())&&(!ss.eof()))
        return FALSE;
      ss>>c; if(c!='-') return FALSE;
      ss>>month;
      if((!ss.good())&&(!ss.eof()))
        return FALSE;
      ss>>c; if(c!='-') return FALSE;
      ss>>day;
      if((!ss.good())&&(!ss.eof()))
        return FALSE;
      ss>>c;
      if(c=='T')
      {
        ss>>hour;
        if((!ss.good())&&(!ss.eof()))
          return FALSE;
        ss>>c; if(c!=':') return FALSE; //throw(Exception(": expected"));
        ss>>min;
        if((!ss.good())&&(!ss.eof()))
          return FALSE;
      }
    }

    x = date2time(year, month, day, hour, min, 0);
    return TRUE;
  }
  catch(std::exception &e)
  {
    return FALSE;
  }
}

/***********************************************/

void GnssStationLog2Platform::checkSinexFile(const FileName &fileName, Platform &platform) const
{
  try
  {
    logStatus<<"read SINEX file <"<<fileName<<">"<<Log::endl;
    Sinex sinex;
    readFileSinex(fileName, sinex);
    Platform platformSinex;

    Bool found = FALSE;
    for(std::string &line : sinex.findBlock("SITE/ID")->lines)
    {
      line.resize(std::max(line.size(), UInt(18)), ' ');
      if(String::upperCase(String::trim(line.substr(1, 4))) != platform.markerName || String::trim(line.substr(9, 9)) != platform.markerNumber)
        continue;
      platformSinex.markerName   = platform.markerName;
      platformSinex.markerNumber = platform.markerNumber;
      if(line.size() >= 75)
      {
        Ellipsoid ellipsoid;
        Double lon = (String::toDouble(line.substr(44, 3)) + String::toDouble(line.substr(48, 2))/60 + String::toDouble(line.substr(51, 4))/3600) * DEG2RAD;
        Double lat = String::toDouble(line.substr(56, 3)) * DEG2RAD;
        lat += (lat < 0 ? -1 : 1) * (String::toDouble(line.substr(60, 2))/60 + String::toDouble(line.substr(63, 4))/3600) * DEG2RAD;
        platformSinex.approxPosition = ellipsoid(Angle(lon), Angle(lat), String::toDouble(line.substr(68, 7)));
      }
      found = TRUE;
    }

    if(!found)
    {
      logWarning<<platform.markerName<<"."<<platform.markerNumber<<" not found in SINEX"<<Log::endl;
      return;
    }

    if(platform.approxPosition.r()<6300e3 && platformSinex.approxPosition.r()>6300e3)
    {
      logWarning<<"using approx position from SINEX file"<<Log::endl;
      platform.approxPosition = platformSinex.approxPosition;
    }

    // check receiver
    // --------------
    std::vector<PlatformGnssReceiver> sinexReceivers;
    for(std::string &line : sinex.findBlock("SITE/RECEIVER")->lines)
    {
      line.resize(std::max(line.size(), UInt(80)), ' ');
      if(String::upperCase(String::trim(line.substr(1, 4))) != platformSinex.markerName)
        continue;
      PlatformGnssReceiver receiver;
      receiver.timeStart = Sinex::str2time(line, 16, FALSE);
      receiver.timeEnd   = Sinex::str2time(line, 29, TRUE);
      receiver.name      = String::trim(line.substr(42, 20));
      receiver.serial    = String::trim(line.substr(63, 5));
      receiver.version   = String::trim(line.substr(69, 11));
      sinexReceivers.emplace_back(receiver);
    }

    for(const auto &equipment : platform.equipments)
    {
      auto recv = std::dynamic_pointer_cast<PlatformGnssReceiver>(equipment);
      if(!recv)
        continue;

      auto sinexRecv = std::lower_bound(sinexReceivers.begin(), sinexReceivers.end(), recv->timeStart,
                                  [](const auto &x, const Time &time) {return x.timeEnd <= time;});
      if((sinexRecv == sinexReceivers.end()) || (sinexRecv->timeStart > recv->timeStart))
      {
        logWarning<<platform.markerName<<"."<<platform.markerNumber<<": receiver sinex check failed (no receiver found in sinex file at time: "
                  <<recv->timeStart.dateTimeStr()<<")"<<Log::endl;
        continue;
      }

      // temporarily combine sinex file entries if platform log spans multiple entries
      for(auto sinexRecv2=sinexRecv+1; (sinexRecv2 != sinexReceivers.end()) && (sinexRecv2->timeStart < recv->timeEnd); sinexRecv2++)
        if((sinexRecv->name    == sinexRecv2->name) &&
           (sinexRecv->serial  == sinexRecv2->serial) &&
           (sinexRecv->version == sinexRecv2->version))
          sinexRecv->timeEnd = sinexRecv2->timeEnd;

      // test receiver name
      if(recv->name != sinexRecv->name)
        logWarning<<platform.markerName<<"."<<platform.markerNumber<<": receiver name sinex check failed (log: "
                  <<recv->name<<", sinex: "<<sinexRecv->name<<")"<<Log::endl;

      // test receiver serial
      if(recv->serial.substr(0,sinexRecv->serial.length()) != sinexRecv->serial)
        logWarning<<platform.markerName<<"."<<platform.markerNumber<<": receiver serial sinex check failed (log: "
                  <<recv->serial<<", sinex: "<<sinexRecv->serial<<")"<<Log::endl;

      // test receiver version
      if(recv->version.substr(0,sinexRecv->version.length()) != sinexRecv->version)
        logWarning<<platform.markerName<<"."<<platform.markerNumber<<": receiver version sinex check failed (log: "
                  <<recv->version<<", sinex: "<<sinexRecv->version<<")"<<Log::endl;

      // test receiver timeStart
      if(std::fabs((recv->timeStart - sinexRecv->timeStart).seconds()) > 1.1)
        logWarning<<platform.markerName<<"."<<platform.markerNumber<<": receiver timeStart sinex check failed (log: "
                  <<recv->timeStart.dateTimeStr()<<", sinex: "<<sinexRecv->timeStart.dateTimeStr() <<")"<<Log::endl;

      // test receiver timeEnd
      if(std::fabs((recv->timeEnd - sinexRecv->timeEnd).seconds()) > 1.1)
        logWarning<<platform.markerName<<"."<<platform.markerNumber<<": receiver timeEnd sinex check failed (log: "
                  <<recv->timeEnd.dateTimeStr()<<", sinex: "<<sinexRecv->timeEnd.dateTimeStr() <<")"<<Log::endl;
    }

    // check antenna
    // -------------
    std::vector<PlatformGnssAntenna> sinexAntennas;
    for(std::string &line : sinex.findBlock("SITE/ANTENNA")->lines)
    {
      line.resize(std::max(line.size(), UInt(68)), ' ');
      if(String::upperCase(String::trim(line.substr(1, 4))) != platformSinex.markerName)
        continue;
      PlatformGnssAntenna antenna;
      antenna.timeStart = Sinex::str2time(line, 16, FALSE);
      antenna.timeEnd   = Sinex::str2time(line, 29, TRUE);
      antenna.name      = String::trim(line.substr(42, 15));
      antenna.radome    = String::trim(line.substr(58, 4));
      antenna.serial    = String::trim(line.substr(63, 5));
      if(antenna.radome == "none" || antenna.radome == "NONE")
        antenna.radome = "";
      if(antenna.serial == "n/a" || antenna.serial == "N/A")
        antenna.serial = "";
      sinexAntennas.emplace_back(antenna);
    }

    for(std::string &line : sinex.findBlock("SITE/ECCENTRICITY")->lines)
    {
      line.resize(std::max(line.size(), UInt(72)), ' ');
      if(String::upperCase(String::trim(line.substr(1, 4))) != platformSinex.markerName)
        continue;
      Time timeStart     = Sinex::str2time(line, 16, FALSE);
      Time timeEnd       = Sinex::str2time(line, 29, TRUE);
      std::string refSys = String::trim(line.substr(42, 3));
      Vector3d eccentricity(String::toDouble(line.substr(55, 8)), String::toDouble(line.substr(64, 8)), String::toDouble(line.substr(46, 8)));

      if(refSys != "UNE")
      {
        logWarning<<platformSinex.markerName<<"."<<platformSinex.markerNumber<<": unknown eccentricity reference system '"<<refSys<<"' in sinex file for time period "
                  <<timeStart.dateTimeStr()<<" to "<<timeEnd.dateTimeStr() <<Log::endl;
        continue;
      }

      Bool found = FALSE;
      for(UInt i=0; i<sinexAntennas.size(); i++)
        if((sinexAntennas.at(i).timeStart < timeEnd) && (sinexAntennas.at(i).timeEnd > timeStart))
        {
          found = TRUE;
          sinexAntennas.at(i).position = eccentricity;
          // split antenna that fits eccentricity time period if eccentricity ends before antenna
          if(sinexAntennas.at(i).timeEnd > timeEnd)
          {
            sinexAntennas.insert(sinexAntennas.begin()+i, sinexAntennas.at(i));
            sinexAntennas.at(i).timeEnd     = timeEnd;
            sinexAntennas.at(i+1).timeStart = timeEnd;
          }
        }

      if(!found)
        logWarning<<platformSinex.markerName<<"."<<platformSinex.markerNumber<<": no antenna found in sinex file to match eccentricity for time period "
                  <<timeStart.dateTimeStr()<<" to "<<timeEnd.dateTimeStr() <<Log::endl;
    }

    // ===========================================

    for(const auto &equipment : platform.equipments)
    {
      auto ant = std::dynamic_pointer_cast<PlatformGnssAntenna>(equipment);
      if(!ant)
        continue;

      auto sinexAnt = std::lower_bound(sinexAntennas.begin(), sinexAntennas.end(), ant->timeStart,
                                  [](const auto &x, const Time &time) {return x.timeEnd <= time;});
      if((sinexAnt == sinexAntennas.end()) || (sinexAnt->timeStart > ant->timeStart))
      {
        logWarning<<platform.markerName<<"."<<platform.markerNumber<<": antenna sinex check failed (no antenna found in sinex file at time: "
                  <<ant->timeStart.dateTimeStr()<<")"<<Log::endl;
        continue;
      }

      // temporarily combine sinex file entries if platform log spans multiple entries
      for(auto sinexAnt2=sinexAnt+1; (sinexAnt2 != sinexAntennas.end()) && (sinexAnt2->timeStart < ant->timeEnd); sinexAnt2++)
        if((sinexAnt->name   == sinexAnt2->name) &&
           (sinexAnt->serial == sinexAnt2->serial) &&
           (sinexAnt->radome == sinexAnt2->radome))
          sinexAnt->timeEnd = sinexAnt2->timeEnd;

      // test antenna name
      if(ant->name != sinexAnt->name)
        logWarning<<platform.markerName<<"."<<platform.markerNumber<<": antenna name sinex check failed (log: "
                  <<ant->name<<", sinex: "<<sinexAnt->name<<")"<<Log::endl;

      // test antenna serial
      if(ant->serial.substr(0,sinexAnt->serial.length()) != sinexAnt->serial)
        logWarning<<platform.markerName<<"."<<platform.markerNumber<<": antenna serial sinex check failed (log: "
                  <<ant->serial<<", sinex: "<<sinexAnt->serial<<")"<<Log::endl;

      // test antenna radome
      if(ant->radome.substr(0,sinexAnt->radome.length()) != sinexAnt->radome)
        logWarning<<platform.markerName<<"."<<platform.markerNumber<<": antenna radome sinex check failed (log: "
                  <<ant->radome<<", sinex: "<<sinexAnt->radome<<")"<<Log::endl;

      // test antenna eccentricity
      if((ant->position-sinexAnt->position).r() > 0.0001)
        logWarning<<platform.markerName<<"."<<platform.markerNumber<<": antenna eccentricity sinex check failed (log: [n,e,u]=["
                  <<ant->position.x()<<","<<ant->position.y()<<","<<ant->position.z()
                  <<"], sinex: [n,e,u]=["<<sinexAnt->position.x()<<","<<sinexAnt->position.y()<<","
                  <<sinexAnt->position.z()<<"])"<<Log::endl;

      // test antenna timeStart
      if(std::fabs((ant->timeStart - sinexAnt->timeStart).seconds()) > 1.1)
        logWarning<<platform.markerName<<"."<<platform.markerNumber<<": antenna timeStart sinex check failed (log: "
                  <<ant->timeStart.dateTimeStr()<<", sinex: "<<sinexAnt->timeStart.dateTimeStr() <<")"<<Log::endl;

      // test antenna timeEnd
      if(std::fabs((ant->timeEnd - sinexAnt->timeEnd).seconds()) > 1.1)
        logWarning<<platform.markerName<<"."<<platform.markerNumber<<": antenna timeEnd sinex check failed (log: "
                  <<ant->timeEnd.dateTimeStr()<<", sinex: "<<sinexAnt->timeEnd.dateTimeStr() <<")"<<Log::endl;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
