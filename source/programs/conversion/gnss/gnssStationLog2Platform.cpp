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
Cross-checking station log data with a \href{https://www.iers.org/IERS/EN/Organization/AnalysisCoordinator/SinexFormat/sinex.html}{SINEX file} is
possible with \program{CheckStationsPlatformsWithSinex}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "inputOutput/file.h"
#include "files/filePlatform.h"

/***** CLASS ***********************************/

/** @brief GNSS analysis.
* @ingroup programsConversionGroup */
class GnssStationLog2Platform
{
  Platform readFile(const FileName &fileName);
  void     newAntenna (Platform &platform, PlatformGnssAntenna &antenna) const;
  void     newReceiver(Platform &platform, PlatformGnssReceiver &receiver) const;
  void     newFreqStandard(Platform &platform, PlatformEquipment &freqStandard) const;
  Bool     readDouble(const std::string &line, Double &x) const;
  Bool     readString(const std::string &line, std::string &x) const;
  Bool     readTime  (const std::string &line, Time &x) const;
  Bool     readEffectiveDate(const std::string &line, Time &xStart, Time &xEnd) const;

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
    FileName fileNameStationPlatform, fileNameStationLog, fileNameAntenna;

    renameDeprecatedConfig(config, "outputfileStationInfo", "outputfileStationPlatform", date2time(2023, 1, 4));

    readConfig(config, "outputfileStationPlatform",  fileNameStationPlatform, Config::MUSTSET,  "", "");
    readConfig(config, "inputfileStationLog",        fileNameStationLog,      Config::MUSTSET,  "", "");
    readConfig(config, "inputfileAntennaDefinition", fileNameAntenna,         Config::OPTIONAL, "{groopsDataDir}/gnss/receiverStation/antennaDefinition/igs/igs20/antennaDefinition_igs20.dat", "used to check antennas");
    if(isCreateSchema(config)) return;

    // ==========================

    logStatus<<"read station log <"<<fileNameStationLog<<">"<<Log::endl;
    Platform platform = readFile(fileNameStationLog);

    // ==========================

    // some tests
    // ----------
    if((platform.markerName.size() != 4) && (platform.markerName.size() != 9))
      throw(Exception(platform.markerName+"."+platform.markerNumber+": marker name should have 4 or 9 letters"));
    if(platform.approxPosition.r() < 6300e3)
      throw(Exception(platform.markerName+"."+platform.markerNumber+": No approx. position given"));
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
          logWarning<<platform.markerName<<"."<<platform.markerNumber<<": no antenna definition found for "<<antenna->str()<<Log::endl;
      }
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
    PlatformEquipment    freqStandard;

    InFile file(fileName);

    enum Section {OTHER     = 999,
                  SITE      = 1,
                  LOCATION  = 2,
                  RECEIVER  = 3,
                  ANTENNA   = 4,
                  FREQUENCY = 6};
    Section section = OTHER;

    while(std::getline(file, line))
    {
      if(line.empty())
        continue;
      if(line[line.length()-1] == '\r')
        line.erase(line.length()-1);

      UInt pos = line.find_first_not_of(" \t\n\r");
      if((pos == std::string::npos) || (pos > 30))
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
      else if(line.at(0) == '6')  // Frequency Standard
      {
        section = FREQUENCY;
        newFreqStandard(platform, freqStandard);
        freqStandard = PlatformEquipment();
      }
      else if(line.at(0) == '7')  // Stop with/after Collocation Information
      {
        break;
      }
      else if((line.at(0) != ' ') && (line.at(0) != '\t')) // other section
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
        if(String::contains(line, "Site Identification of the GNSS Monument"))
        {
        }
        else if(String::contains(line, "Four Character ID") || String::contains(line, "Nine Character ID"))
        {
          if(platform.markerName.empty())
          {
            readString(line, platform.markerName);
            std::transform(platform.markerName.begin(), platform.markerName.end(), platform.markerName.begin(), ::toupper);
          }
        }
        else if(String::contains(line, "IERS DOMES Number"))
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
        if(String::contains(line, "Site Location Information"))
        {
        }
        else if(String::contains(line, "City or Town"))
        {
          readString(line, platform.comment);
        }
        else if(String::contains(line, "Country"))
        {
          std::string country;
          readString(line, country);
          platform.comment += ", " + country;
        }
        else if(String::contains(line, "X coordinate"))
        {
          readDouble(line, platform.approxPosition.x());
        }
        else if(String::contains(line, "Y coordinate"))
        {
          readDouble(line, platform.approxPosition.y());
        }
        else if(String::contains(line, "Z coordinate"))
        {
          readDouble(line, platform.approxPosition.z());
        }
        else if(platform.approxPosition.r() < 6300e3 && String::contains(line, "Latitude"))
        {
          readDouble(line, platform.approxPosition.x());
        }
        else if(platform.approxPosition.r() < 6300e3 && String::contains(line, "Longitude"))
        {
          readDouble(line, platform.approxPosition.y());
        }
        else if(platform.approxPosition.r() < 6300e3 && String::contains(line, "Elevation"))
        {
          readDouble(line, platform.approxPosition.z());
        }
/*        else
          logInfo<<"unknown line '"<<line<<"'"<<Log::endl;*/
      } // if(section == LOCATION)

      // =================================================

      if(section == RECEIVER)
      {
        if((String::contains(line, "Additional Information")) || (firstLetter==':'))
        {
//           std::string comment;
//           readString(line, comment);
//           receiver.comment += comment;
        }
        else if(String::contains(line, "GNSS Receiver Information"))
        {
        }
        else if(String::contains(line, "Receiver Type"))
        {
          readString(line, receiver.name);
          if((receiver.name != "ASHTECH Z-XII3T") && String::contains(receiver.name, "ASHTECH") && String::contains(receiver.name, "Z-XII"))
            receiver.name = "ASHTECH Z-XII3";
        }
        else if(String::contains(line, "Satellite System"))
        {
        }
        else if(String::contains(line, "Serial Number"))
        {
          readString(line, receiver.serial);
        }
        else if(String::contains(line, "Firmware Version"))
        {
          readString(line, receiver.version);
        }
        else if(String::contains(line, "Elevation Cutoff"))
        {
        }
        else if(String::contains(line, "Date Installed"))
        {
          readTime(line, receiver.timeStart);
        }
        else if(String::contains(line, "Date Removed"))
        {
          readTime(line, receiver.timeEnd);
        }
        else if(String::contains(line, "Temperature Stabiliz"))
        {
        }
//         else
//           logWarning<<platform.markerName<<"."<<platform.markerNumber<<" unknown line in GNSS Receiver Information section '"<<line<<"'"<<Log::endl;
      } // if(section == RECEIVER)

      // =================================================

      if(section == ANTENNA)
      {
        if((String::contains(line, "Additional Information")) || (firstLetter==':'))
        {
//           std::string comment;
//           readString(line, comment);
//           antenna.comment += comment;
        }
        else if(String::contains(line, "GNSS Antenna Information"))
        {
        }
        else if(String::contains(line, "Antenna Type"))
        {
          std::string text;
          readString(line, text);
          text = String::upperCase(text);
          if(text.size() == 20) // with radome?
          {
            antenna.radome = text.substr(16, 4);
            if(antenna.radome == "NONE")
              antenna.radome.clear();
          }
          // if(text.substr(0, 8) == "TRIMBLE ")
          //   text = text.substr(8);
          if(text == "DORNE MARGOLIN T")
            text = "AOAD/M_T";
          text = String::trim(text.substr(0, 15));
          antenna.name = text;
        }
        else if(String::contains(line, "Radome Serial Number"))
        {
        }
        else if(String::contains(line, "Serial Number"))
        {
          readString(line, antenna.serial);
          antenna.serial = String::upperCase(antenna.serial).substr(0, 5);
          if((antenna.serial=="n/a")||(antenna.serial=="N/A"))
            antenna.serial.clear();
        }
        else if(String::contains(line, "Antenna Reference Point"))
        {
        }
        else if(String::contains(line, "Marker->ARP Up"))
        {
          readDouble(line, antenna.position.z());
        }
        else if(String::contains(line, "Marker->ARP East"))
        {
          readDouble(line, antenna.position.y());
        }
        else if(String::contains(line, "Marker->ARP North"))
        {
          readDouble(line, antenna.position.x());
        }
        else if(String::contains(line, "Alignment from True N") || String::contains(line, "Degree Offset from North") /*old format*/)
        {
          Double deg;
          if(readDouble(line, deg) && (deg!=0))
            antenna.local2antennaFrame = rotaryZ(Angle(deg*DEG2RAD));
        }
        else if(String::contains(line, "Antenna Radome Type"))
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
        else if(String::contains(line, "Antenna Cable Type"))
        {
        }
        else if(String::contains(line, "Antenna Cable Length"))
        {
        }
        else if(String::contains(line, "Date Installed"))
        {
          readTime(line, antenna.timeStart);
        }
        else if(String::contains(line, "Date Removed"))
        {
          readTime(line, antenna.timeEnd);
        }
        else if(String::contains(line, "Antenna Height")) // old format
        {
          readDouble(line, antenna.position.z());
        }
//         else
//           logWarning<<platform.markerName<<"."<<platform.markerNumber<<" unknown line in GNSS Antenna Information section '"<<line<<"'"<<Log::endl;
      } // if(section == ANTENNA)

      // =================================================

      if(section == FREQUENCY)
      {
        if((String::contains(line, "Notes")) || (firstLetter==':'))
        {
        }
        else if(String::contains(line, "Standard Type"))
        {
          std::string text;
          readString(line, text);
          text = String::upperCase(text);
          if(String::contains(text, "H-MASER") || String::contains(text, "HYDROGEN"))
            freqStandard.name = "H-MASER";
        }
        else if(String::contains(line, "Input Frequency"))
        {
        }
        else if(String::contains(line, "Effective Dates"))
        {
          readEffectiveDate(line, freqStandard.timeStart, freqStandard.timeEnd);
        }
      } // if(section == FREQUENCY)

    } // for(;;)

    // if no Cartesian coordinates are given, try to use ellipsoidal coordinates
    if(platform.approxPosition.r() > 0 && platform.approxPosition.r() < 6300e3)
    {
      Ellipsoid ellipsoid;
      platform.approxPosition = ellipsoid(Angle(platform.approxPosition.y()*DEG2RAD), Angle(platform.approxPosition.x()*DEG2RAD), platform.approxPosition.z());
    }

    // insert last receiver/antenna/freqStandard
    // -----------------------------------------
    newReceiver(platform, receiver);
    newAntenna(platform, antenna);
    newFreqStandard(platform, freqStandard);

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
    if(receiver.timeStart == receiver.timeEnd)
      return;
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
    if(antenna.timeStart == antenna.timeEnd)
      return;
    if(antenna.timeStart >= antenna.timeEnd)
      logWarning<<platform.markerName<<"."<<platform.markerNumber<<": antenna time error (new antenna start: "
                <<antenna.timeStart.dateTimeStr()<<", new antenna end: "<<antenna.timeEnd.dateTimeStr()<<")"<<Log::endl;

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
       ((lastAntenna->position-antenna.position).r() < 0.0001) &&
      quadsum(lastAntenna->local2antennaFrame.matrix() - antenna.local2antennaFrame.matrix()) < 0.0001)
    {
      lastAntenna->timeEnd = antenna.timeEnd;
      return;
    }

    // test times
    if(lastAntenna->timeEnd == date2time(2500,1,1))
      lastAntenna->timeEnd = antenna.timeStart;

    if(lastAntenna->timeEnd > antenna.timeStart)
      logWarning<<platform.markerName<<"."<<platform.markerNumber<<" antenna time error (previous antenna end: "
                <<lastAntenna->timeEnd.dateTimeStr()<<", new antenna start: "<<antenna.timeStart.dateTimeStr()<<")"<<Log::endl;

    platform.equipments.push_back(std::make_shared<PlatformGnssAntenna>(antenna));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssStationLog2Platform::newFreqStandard(Platform &platform, PlatformEquipment &freqStandard) const
{
  try
  {
    if(freqStandard.name.empty())
      return;

    if(freqStandard.timeEnd == Time())
      freqStandard.timeEnd = date2time(2500,1,1);
    else
      freqStandard.timeEnd += mjd2time(1);
    if(freqStandard.timeStart == freqStandard.timeEnd)
      return;
    if(freqStandard.timeStart >= freqStandard.timeEnd)
      logWarning<<platform.markerName<<"."<<platform.markerNumber<<": freqStandard time error (new freqStandard start: "
                <<freqStandard.timeStart.dateTimeStr()<<", new freqStandard end: "<<freqStandard.timeEnd.dateTimeStr()<<")"<<Log::endl;

    // find last freqStandard
    auto iter = std::find_if(platform.equipments.rbegin(), platform.equipments.rend(),
                             [](const auto x) {return x->getType() == PlatformEquipment::OTHER;});
    if(iter == platform.equipments.rend())
    {
      platform.equipments.push_back(std::make_shared<PlatformEquipment>(freqStandard));
      return;
    }
    auto lastFreqStandard = std::dynamic_pointer_cast<PlatformEquipment>(*iter);

    if((lastFreqStandard->name == freqStandard.name) && (lastFreqStandard->serial == freqStandard.serial))
    {
      lastFreqStandard->timeEnd = freqStandard.timeEnd;
      return;
    }

    // test times
    if(lastFreqStandard->timeEnd == date2time(2500,1,1))
      lastFreqStandard->timeEnd = freqStandard.timeStart;

    if(lastFreqStandard->timeEnd > freqStandard.timeStart)
      logWarning<<platform.markerName<<"."<<platform.markerNumber<<" freqStandard time error (previous freqStandard end: "
                <<lastFreqStandard->timeEnd.dateTimeStr()<<", new freqStandard start: "<<freqStandard.timeStart.dateTimeStr()<<")"<<Log::endl;

    platform.equipments.push_back(std::make_shared<PlatformEquipment>(freqStandard));
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
    std::string text = String::trim(line.substr(pos+1));
    if(text.empty() || String::startsWith(text, "("))
      return FALSE;
    x = text;
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
    std::string text = String::trim(line.substr(pos+1));
    if(text.empty() || String::startsWith(text, "("))
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
      ss>>hour>>c>>min;  // ignore problems when reading the clock
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
        ss>>hour>>c>>min; // ignore problems when reading the clock
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

Bool GnssStationLog2Platform::readEffectiveDate(const std::string &line, Time &xStart, Time &xEnd) const
{
  try
  {
    UInt pos = line.find(":");
    if(pos == std::string::npos)
      throw(Exception("':' expected"));
    std::string text = String::trim(line.substr(pos+1));
    if(text.empty() || String::startsWith(text, "("))
      return FALSE;
    std::stringstream ss(text);

    UInt year, month, day;
    char c;
    ss>>year;
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
    xStart = date2time(year, month, day);

    ss>>c; if(c!='/') return TRUE;
    ss>>year>>c>>month>>c>>day;
    if((!ss.good())&&(!ss.eof()))
      return TRUE;
    xEnd = date2time(year, month, day);

    return TRUE;
  }
  catch(std::exception &e)
  {
    return FALSE;
  }
}

/***********************************************/
