/***********************************************/
/**
* @file gnssStationLog2StationInfo.cpp
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
Converts \href{https://files.igs.org/pub/station/general/blank.log}{IGS station log format} to \configFile{outputfileStationInfo}{gnssStationInfo}.

If \configFile{inputfileAntennaDefinition}{gnssAntennaDefinition} is provided, station log data is cross-checked with the given antenna definitions.
Cross-checking station log data with a \href{https://www.iers.org/IERS/EN/Organization/AnalysisCoordinator/SinexFormat/sinex.html}{SINEX file} is also
possible by providing \config{inputfileSinex}. Any failed checks result in warnings in the output log.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "inputOutput/file.h"
#include "inputOutput/fileSinex.h"
#include "files/fileGnssStationInfo.h"

/***** CLASS ***********************************/

/** @brief GNSS analysis.
* @ingroup programsConversionGroup */
class GnssStationLog2StationInfo
{
  GnssStationInfo readFile(const FileName &fileName);
  GnssStationInfo readSinexFile(const FileName &fileName, const GnssStationInfo &station) const;
  void            readSinexStation(std::string &line, const GnssStationInfo& station, GnssStationInfo &sinexStationInfo) const;
  void            readSinexReceiver(std::string &line, GnssStationInfo &sinexStationInfo) const;
  void            readSinexAntenna(std::string &line, GnssStationInfo &sinexStationInfo) const;
  void            readSinexEccentricity(std::string &line, GnssStationInfo &sinexStationInfo) const;
  void            newAntenna (GnssStationInfo &station, GnssAntennaInfo &antenna) const;
  void            newReceiver(GnssStationInfo &station, GnssReceiverInfo &receiver) const;
  Bool            readDouble(const std::string &line, Double &x) const;
  Bool            readString(const std::string &line, std::string &x) const;
  Bool            readTime  (const std::string &line, Time &x) const;

public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GnssStationLog2StationInfo, SINGLEPROCESS, "GNSS analysis", Conversion, Gnss)

/***********************************************/

void GnssStationLog2StationInfo::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameStationLog, fileNameAntenna, fileNameStationInfo, fileNameSinex;

    readConfig(config, "outputfileStationInfo",      fileNameStationInfo, Config::MUSTSET,   "", "");
    readConfig(config, "inputfileStationLog",        fileNameStationLog,  Config::MUSTSET,   "", "");
    readConfig(config, "inputfileAntennaDefinition", fileNameAntenna,     Config::OPTIONAL,  "{groopsDataDir}/gnss/receiverStation/antennaDefinition/igs/igs20/antennaDefinition_igs20.dat", "used to check antennas");
    readConfig(config, "inputfileSinex",             fileNameSinex,       Config::OPTIONAL,  "{groopsDataDir}/gnss/receiverStation/igs_with_former.snx", "used to cross-check station log with SINEX file");
    if(isCreateSchema(config)) return;

    // ==========================

    logStatus<<"read station log <"<<fileNameStationLog<<">"<<Log::endl;
    GnssStationInfo station = readFile(fileNameStationLog);

    // ==========================

    // some tests
    // ----------
    if(station.markerName.size()!=4)
    {
      logWarning<<station.markerName<<"."<<station.markerNumber<<": marker name should have 4 letters"<<Log::endl;
      throw(Exception("marker name should have 4 letters"));
    }
    if(station.receiver.size()==0)
    {
      logWarning<<station.markerName<<"."<<station.markerNumber<<": No receiver given"<<Log::endl;
      throw(Exception("No receiver given"));
    }
    if(station.antenna.size()==0)
    {
      logWarning<<station.markerName<<"."<<station.markerNumber<<": No antenna given"<<Log::endl;
      throw(Exception("No antenna given"));
    }

    // ==========================

    // test antennas
    // -------------
    if(!fileNameAntenna.empty())
    {
      logStatus<<"read antenna definition file <"<<fileNameAntenna<<">"<<Log::endl;
      std::vector<GnssAntennaDefinitionPtr> antennaList;
      readFileGnssAntennaDefinition(fileNameAntenna, antennaList);

      station.fillAntennaPattern(antennaList);
      // test antennas
      for(UInt idAnt=0; idAnt<station.antenna.size(); idAnt++)
        if(!station.antenna.at(idAnt).antennaDef)
          logWarning<<station.markerName<<"."<<station.markerNumber<<": no antenna definition found for "
                    <<station.antenna.at(idAnt).name+"."<<station.antenna.at(idAnt).serial<<Log::endl;
    }

    // ==========================

    // check station with sinex file
    // -----------------------------
    if(!fileNameSinex.empty())
    {
      logStatus<<"read SINEX file <"<<fileNameSinex<<">"<<Log::endl;
      GnssStationInfo sinexStation = readSinexFile(fileNameSinex, station);

      if(station.approxPosition.r()<6300e3 && sinexStation.approxPosition.r()>6300e3)
      {
        logWarning<<"using approx position from SINEX file"<<Log::endl;
        station.approxPosition = sinexStation.approxPosition;
      }

      if(station.receiver.size() != sinexStation.receiver.size())
        logWarning<<station.markerName<<"."<<station.markerNumber<<": receiver count sinex mismatch: "<<station.receiver.size()<<" != "<<sinexStation.receiver.size()<<Log::endl;

      // receiver check
      for(UInt idx = 0; idx < station.receiver.size(); idx++)
      {
        UInt idxSinex    = sinexStation.findReceiver(station.receiver.at(idx).timeStart);
        UInt idxSinexEnd = sinexStation.findReceiver(station.receiver.at(idx).timeEnd - seconds2time(0.1));

        if(idxSinex == NULLINDEX)
        {
          logWarning<<station.markerName<<"."<<station.markerNumber<<": receiver sinex check failed (no receiver found in sinex file at time: "
                    <<station.receiver.at(idx).timeStart.dateTimeStr()<<")"<<Log::endl;
          continue;
        }

        // temporarily combine sinex file entries if station log spans multiple entries
        if(idxSinexEnd != NULLINDEX && idxSinex < idxSinexEnd)
          for(UInt i = idxSinex+1; i <= idxSinexEnd; i++)
            if(sinexStation.receiver.at(idxSinex).name    == sinexStation.receiver.at(i).name &&
               sinexStation.receiver.at(idxSinex).serial  == sinexStation.receiver.at(i).serial &&
               sinexStation.receiver.at(idxSinex).version == sinexStation.receiver.at(i).version)
              sinexStation.receiver.at(idxSinex).timeEnd = sinexStation.receiver.at(i).timeEnd;

        // test receiver name
        if(station.receiver.at(idx).name != sinexStation.receiver.at(idxSinex).name)
          logWarning<<station.markerName<<"."<<station.markerNumber<<": receiver name sinex check failed (log: "
                    <<station.receiver.at(idx).name<<", sinex: "<<sinexStation.receiver.at(idxSinex).name<<")"<<Log::endl;

        // test receiver serial
        if(station.receiver.at(idx).serial.substr(0,sinexStation.receiver.at(idxSinex).serial.length()) != sinexStation.receiver.at(idxSinex).serial)
          logWarning<<station.markerName<<"."<<station.markerNumber<<": receiver serial sinex check failed (log: "
                    <<station.receiver.at(idx).serial<<", sinex: "<<sinexStation.receiver.at(idxSinex).serial<<")"<<Log::endl;

        // test receiver version
        if(station.receiver.at(idx).version.substr(0,sinexStation.receiver.at(idxSinex).version.length()) != sinexStation.receiver.at(idxSinex).version)
          logWarning<<station.markerName<<"."<<station.markerNumber<<": receiver version sinex check failed (log: "
                    <<station.receiver.at(idx).version<<", sinex: "<<sinexStation.receiver.at(idxSinex).version<<")"<<Log::endl;

        // test receiver timeStart
        if(station.receiver.at(idx).timeStart != sinexStation.receiver.at(idxSinex).timeStart &&
           fabs((station.receiver.at(idx).timeStart - sinexStation.receiver.at(idxSinex).timeStart).seconds()) > 1.1)
          logWarning<<station.markerName<<"."<<station.markerNumber<<": receiver timeStart sinex check failed (log: "
                    <<station.receiver.at(idx).timeStart.dateTimeStr()<<", sinex: "<<sinexStation.receiver.at(idxSinex).timeStart.dateTimeStr() <<")"<<Log::endl;

        // test receiver timeEnd
        if(station.receiver.at(idx).timeEnd != sinexStation.receiver.at(idxSinex).timeEnd &&
           fabs((station.receiver.at(idx).timeEnd - sinexStation.receiver.at(idxSinex).timeEnd).seconds()) > 1.1)
          logWarning<<station.markerName<<"."<<station.markerNumber<<": receiver timeEnd sinex check failed (log: "
                    <<station.receiver.at(idx).timeEnd.dateTimeStr()<<", sinex: "<<sinexStation.receiver.at(idxSinex).timeEnd.dateTimeStr() <<")"<<Log::endl;
      }

      if(station.antenna.size() != sinexStation.antenna.size())
        logWarning<<station.markerName<<"."<<station.markerNumber<<": antenna count sinex mismatch: "<<station.receiver.size()<<" != "<<sinexStation.receiver.size()<<Log::endl;

      // antenna check
      for(UInt idx = 0; idx < station.antenna.size(); idx++)
      {
        UInt idxSinex    = sinexStation.findAntenna(station.antenna.at(idx).timeStart);
        UInt idxSinexEnd = sinexStation.findAntenna(station.antenna.at(idx).timeEnd - seconds2time(0.1));

        if(idxSinex == NULLINDEX)
        {
          logWarning<<station.markerName<<"."<<station.markerNumber<<": antenna sinex check failed (no antenna found in sinex file at time: "
                    <<station.antenna.at(idx).timeStart.dateTimeStr()<<")"<<Log::endl;
          continue;
        }

        // temporarily combine sinex file entries if station log spans multiple entries
        if(idxSinexEnd != NULLINDEX && idxSinex < idxSinexEnd)
          for(UInt i = idxSinex+1; i <= idxSinexEnd; i++)
            if(sinexStation.antenna.at(idxSinex).name    == sinexStation.antenna.at(i).name &&
               sinexStation.antenna.at(idxSinex).serial  == sinexStation.antenna.at(i).serial &&
               sinexStation.antenna.at(idxSinex).radome == sinexStation.antenna.at(i).radome)
              sinexStation.antenna.at(idxSinex).timeEnd = sinexStation.antenna.at(i).timeEnd;

        // test antenna name
        if(station.antenna.at(idx).name != sinexStation.antenna.at(idxSinex).name)
          logWarning<<station.markerName<<"."<<station.markerNumber<<": antenna name sinex check failed (log: "
                    <<station.antenna.at(idx).name<<", sinex: "<<sinexStation.antenna.at(idxSinex).name<<")"<<Log::endl;

        // test antenna serial
        if(station.antenna.at(idx).serial.substr(0,sinexStation.antenna.at(idxSinex).serial.length()) != sinexStation.antenna.at(idxSinex).serial)
          logWarning<<station.markerName<<"."<<station.markerNumber<<": antenna serial sinex check failed (log: "
                    <<station.antenna.at(idx).serial<<", sinex: "<<sinexStation.antenna.at(idxSinex).serial<<")"<<Log::endl;

        // test antenna radome
        if(station.antenna.at(idx).radome.substr(0,sinexStation.antenna.at(idxSinex).radome.length()) != sinexStation.antenna.at(idxSinex).radome)
          logWarning<<station.markerName<<"."<<station.markerNumber<<": antenna radome sinex check failed (log: "
                    <<station.antenna.at(idx).radome<<", sinex: "<<sinexStation.antenna.at(idxSinex).radome<<")"<<Log::endl;

        // test antenna eccentricity
        if(round(station.antenna.at(idx).position.x()*1e4) != round(sinexStation.antenna.at(idxSinex).position.x()*1e4) ||
           round(station.antenna.at(idx).position.y()*1e4) != round(sinexStation.antenna.at(idxSinex).position.y()*1e4) ||
           round(station.antenna.at(idx).position.z()*1e4) != round(sinexStation.antenna.at(idxSinex).position.z()*1e4))
          logWarning<<station.markerName<<"."<<station.markerNumber<<": antenna eccentricity sinex check failed (log: [n,e,u]=["
                    <<station.antenna.at(idx).position.x()<<","<<station.antenna.at(idx).position.y()<<","<<station.antenna.at(idx).position.z()
                    <<"], sinex: [n,e,u]=["<<sinexStation.antenna.at(idxSinex).position.x()<<","<<sinexStation.antenna.at(idxSinex).position.y()<<","
                    <<sinexStation.antenna.at(idxSinex).position.z()<<"])"<<Log::endl;

        // test antenna timeStart
        if(station.antenna.at(idx).timeStart != sinexStation.antenna.at(idxSinex).timeStart &&
           fabs((station.antenna.at(idx).timeStart - sinexStation.antenna.at(idxSinex).timeStart).seconds()) > 1.1)
          logWarning<<station.markerName<<"."<<station.markerNumber<<": antenna timeStart sinex check failed (log: "
                    <<station.antenna.at(idx).timeStart.dateTimeStr()<<", sinex: "<<sinexStation.antenna.at(idxSinex).timeStart.dateTimeStr() <<")"<<Log::endl;

        // test antenna timeEnd
        if(station.antenna.at(idx).timeEnd != sinexStation.antenna.at(idxSinex).timeEnd &&
           fabs((station.antenna.at(idx).timeEnd - sinexStation.antenna.at(idxSinex).timeEnd).seconds()) > 1.1)
          logWarning<<station.markerName<<"."<<station.markerNumber<<": antenna timeEnd sinex check failed (log: "
                    <<station.antenna.at(idx).timeEnd.dateTimeStr()<<", sinex: "<<sinexStation.antenna.at(idxSinex).timeEnd.dateTimeStr() <<")"<<Log::endl;
      }
    }

    // test if there is still no approx position (neither in station log nor SINEX file)
    if(station.approxPosition.r()<6300e3)
    {
      logWarning<<station.markerName<<"."<<station.markerNumber<<": No approx. position given"<<Log::endl;
      throw(Exception("No approx. position given"));
    }

    // ==========================

    logStatus<<"write station info <"<<fileNameStationInfo<<">"<<Log::endl;
    writeFileGnssStationInfo(fileNameStationInfo, station);

    // test station
    try
    {
      readFileGnssStationInfo(fileNameStationInfo, station);
    }
    catch(std::exception &/*e*/)
    {
      logWarning<<"not readable: '"<<fileNameStationInfo<<"'"<<Log::endl;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

GnssStationInfo GnssStationLog2StationInfo::readFile(const FileName &fileName)
{
  std::string line;

  try
  {
    GnssStationInfo  station;
    GnssReceiverInfo receiver;
    GnssAntennaInfo  antenna;

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

        newReceiver(station, receiver);
        receiver = GnssReceiverInfo();
      }
      else if(line.at(0) == '4')  // GNSS Antenna Information
      {
        section = ANTENNA;

        newAntenna(station, antenna);
        antenna = GnssAntennaInfo();
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
          if(station.markerName.empty())
          {
            readString(line, station.markerName);
            std::transform(station.markerName.begin(), station.markerName.end(), station.markerName.begin(), ::toupper);
            if(station.markerName.length() > 4)
            {
              logWarning<<station.markerName<<": marker name has more than 4 letters, shortened to: "<<station.markerName.substr(0,4)<<Log::endl;
              station.markerName = station.markerName.substr(0,4);
            }
          }
        }
        else if(line.find("IERS DOMES Number") != std::string::npos)
        {
          if(station.markerNumber.empty())
            readString(line, station.markerNumber);
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
          readString(line, station.comment);
        }
        else if(line.find("Country") != std::string::npos)
        {
          std::string country;
          readString(line, country);
          station.comment += ", " + country;
        }
        else if(line.find("X coordinate") != std::string::npos)
        {
          readDouble(line, station.approxPosition.x());
        }
        else if(line.find("Y coordinate") != std::string::npos)
        {
          readDouble(line, station.approxPosition.y());
        }
        else if(line.find("Z coordinate") != std::string::npos)
        {
          readDouble(line, station.approxPosition.z());
        }
        else if(station.approxPosition.r() < 6300e3 && line.find("Latitude") != std::string::npos)
        {
          readDouble(line, station.approxPosition.x());
        }
        else if(station.approxPosition.r() < 6300e3 && line.find("Longitude") != std::string::npos)
        {
          readDouble(line, station.approxPosition.y());
        }
        else if(station.approxPosition.r() < 6300e3 && line.find("Elevation") != std::string::npos)
        {
          readDouble(line, station.approxPosition.z());
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
//           logWarning<<station.markerName<<"."<<station.markerNumber<<" unknown line in GNSS Receiver Information section '"<<line<<"'"<<Log::endl;
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
//           logWarning<<station.markerName<<"."<<station.markerNumber<<" unknown line in GNSS Antenna Information section '"<<line<<"'"<<Log::endl;
      } // if(section == ANTENNA)

    } // for(;;)

    // if no Cartesian coordinates are given, try to use ellipsoidal coordinates
    if(station.approxPosition.r() > 0 && station.approxPosition.r() < 6300e3)
    {
      Ellipsoid ellipsoid;
      station.approxPosition = ellipsoid(Angle(station.approxPosition.y()*DEG2RAD), Angle(station.approxPosition.x()*DEG2RAD), station.approxPosition.z());
    }

    // insert last receiver
    // --------------------
    newReceiver(station, receiver);

    // insert last antenna
    // -------------------
    newAntenna(station, antenna);

    return station;
  }
  catch(std::exception &e)
  {
    logError<<"error: '"<<line<<"'"<<Log::endl;
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssStationLog2StationInfo::newReceiver(GnssStationInfo &station, GnssReceiverInfo &receiver) const
{
  try
  {
    if(receiver.name.empty() || receiver.name == "no receiver")
      return;

    if(receiver.timeStart == Time())
      logWarning<<station.markerName<<"."<<station.markerNumber<<" receiver timeStart not given"<<Log::endl;

    if(receiver.timeEnd == Time())
      receiver.timeEnd = date2time(2500,1,1);
    if(receiver.timeStart == receiver.timeEnd)
      return;
    if(receiver.timeStart >= receiver.timeEnd)
      logWarning<<station.markerName<<"."<<station.markerNumber<<": receiver time error (new receiver start: "
                <<receiver.timeStart.dateTimeStr()<<", new receiver end: "<<receiver.timeEnd.dateTimeStr()<<")"<<Log::endl;

    if(station.receiver.size()==0)
    {
      station.receiver.push_back(receiver);
      return;
    }

    if((station.receiver.back().name    == receiver.name)   &&
       (station.receiver.back().serial  == receiver.serial) &&
       (station.receiver.back().version == receiver.version))
    {
//       logWarning<<station.markerName<<"."<<station.markerNumber<<" receiver not changed"<<Log::endl;
      station.receiver.back().timeEnd = receiver.timeEnd;
      return;
    }

    // test times
    if(station.receiver.back().timeEnd == date2time(2500,1,1))
      station.receiver.back().timeEnd = receiver.timeStart - seconds2time(1.);

    if(station.receiver.back().timeEnd > receiver.timeStart)
      logWarning<<station.markerName<<"."<<station.markerNumber<<" receiver time error (previous receiver end: "
               <<station.receiver.back().timeEnd.dateTimeStr()<<", new receiver start: "<<receiver.timeStart.dateTimeStr()<<")"<<Log::endl;

    station.receiver.push_back(receiver);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssStationLog2StationInfo::newAntenna(GnssStationInfo &station, GnssAntennaInfo &antenna) const
{
  try
  {
    if(antenna.name.empty() || antenna.name == "no antenna")
      return;

    if(antenna.timeStart == Time())
      logWarning<<station.markerName<<"."<<station.markerNumber<<" antenna timeStart not given"<<Log::endl;

    if(antenna.timeEnd == Time())
      antenna.timeEnd = date2time(2500,1,1);
    if(antenna.timeStart == antenna.timeEnd)
      return;
    if(antenna.timeStart >= antenna.timeEnd)
      logWarning<<station.markerName<<"."<<station.markerNumber<<": antenna time error (new antenna start: "
                <<antenna.timeStart.dateTimeStr()<<", new antenna end: "<<antenna.timeEnd.dateTimeStr()<<")"<<Log::endl;

    if(station.antenna.size()==0)
    {
      station.antenna.push_back(antenna);
      return;
    }

    if((station.antenna.back().name    == antenna.name)   &&
       (station.antenna.back().serial  == antenna.serial) &&
       (station.antenna.back().radome  == antenna.radome) &&
       ((station.antenna.back().position-antenna.position).r()<0.0001) &&
       quadsum(station.antenna.back().local2antennaFrame.matrix() - antenna.local2antennaFrame.matrix())<0.0001)
    {
//       logWarning<<station.markerName<<"."<<station.markerNumber<<" antenna not changed"<<Log::endl;
      station.antenna.back().timeEnd = antenna.timeEnd;
      return;
    }

    // test times
    if(station.antenna.back().timeEnd == date2time(2500,1,1))
      station.antenna.back().timeEnd = antenna.timeStart - seconds2time(1.);

    if(station.antenna.back().timeEnd > antenna.timeStart)
      logWarning<<station.markerName<<"."<<station.markerNumber<<" antenna time error (previous antenna end: "
                <<station.antenna.back().timeEnd.dateTimeStr()<<", new antenna start: "<<antenna.timeStart.dateTimeStr()<<")"<<Log::endl;


    station.antenna.push_back(antenna);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool GnssStationLog2StationInfo::readDouble(const std::string &line, Double &x) const
{
  try
  {
    UInt pos = line.find(":");
    if(pos == std::string::npos)
      throw(Exception("':' expected"));
    std::string text = line.substr(pos+1);
//     if(text.find("(") != std::string::npos)
//       return FALSE;
    std::stringstream ss(text);
    Double f = 0;
    ss>>f;
    if((!ss.good())&&(!ss.eof()))
      return FALSE;
    x = f;
// logInfo<<"read double: '"<<x<<"'"<<Log::endl;
    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool GnssStationLog2StationInfo::readString(const std::string &line, std::string &x) const
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
// logInfo<<"read string: '"<<x<<"'"<<Log::endl;
    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool GnssStationLog2StationInfo::readTime(const std::string &line, Time &x) const
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
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

GnssStationInfo GnssStationLog2StationInfo::readSinexFile(const FileName &fileName, const GnssStationInfo &station) const
{
  Sinex sinex;
  readFileSinex(fileName, sinex);
  GnssStationInfo sinexStationInfo;

  for(std::string &line : sinex.findBlock("SITE/ID")->lines)
    readSinexStation(line, station, sinexStationInfo);

  for(std::string &line : sinex.findBlock("SITE/RECEIVER")->lines)
    readSinexReceiver(line, sinexStationInfo);

  for(std::string &line : sinex.findBlock("SITE/ANTENNA")->lines)
    readSinexAntenna(line, sinexStationInfo);

  for(std::string &line : sinex.findBlock("SITE/ECCENTRICITY")->lines)
    readSinexEccentricity(line, sinexStationInfo);

  return sinexStationInfo;
}

/***********************************************/

// +SITE/ID
void GnssStationLog2StationInfo::readSinexStation(std::string &line, const GnssStationInfo &station, GnssStationInfo &sinexStationInfo) const
{
  try
  {
    line.resize(std::max(line.size(), UInt(18)), ' ');
    std::string site = String::trim(line.substr(1, 4));
    std::string id   = String::trim(line.substr(9, 9));
    std::transform(site.begin(), site.end(), site.begin(), ::toupper);

    if(site == station.markerName && id == station.markerNumber)
    {
      sinexStationInfo.markerName   = station.markerName;
      sinexStationInfo.markerNumber = station.markerNumber;
      if(line.size() >= 75)
      {
        const Double longitude   = String::toDouble(line.substr(44, 3)) + String::toDouble(line.substr(48, 2))/60 + String::toDouble(line.substr(51, 4))/3600;
        const Double latitude    = String::toDouble(line.substr(56, 3)) + (String::startsWith(String::trim(line.substr(56, 3)), "-") ? -1 : 1) * (String::toDouble(line.substr(60, 2))/60 + String::toDouble(line.substr(63, 4))/3600);
        const Double height      = String::toDouble(line.substr(68, 7));
        sinexStationInfo.approxPosition = Ellipsoid()(Angle(DEG2RAD*longitude), Angle(DEG2RAD*latitude), height);
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

// +SITE/RECEIVER
void GnssStationLog2StationInfo::readSinexReceiver(std::string &line, GnssStationInfo &station) const
{
  try
  {
    line.resize(std::max(line.size(), UInt(80)), ' ');
    std::string site = String::trim(line.substr(1, 4));
    std::transform(site.begin(), site.end(), site.begin(), ::toupper);

    if(site == station.markerName)
    {
      GnssReceiverInfo receiver;
      receiver.timeStart = Sinex::str2time(line, 16, FALSE);
      receiver.timeEnd   = Sinex::str2time(line, 29, TRUE);
      receiver.name      = String::trim(line.substr(42, 20));
      receiver.serial    = String::trim(line.substr(63, 5));
      receiver.version   = String::trim(line.substr(69, 11));

      station.receiver.push_back(receiver);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

// +SITE/ANTENNA
void GnssStationLog2StationInfo::readSinexAntenna(std::string &line, GnssStationInfo &station) const
{
  try
  {
    line.resize(std::max(line.size(), UInt(68)), ' ');
    std::string site = String::trim(line.substr(1, 4));
    std::transform(site.begin(), site.end(), site.begin(), ::toupper);

    if(site == station.markerName)
    {
      GnssAntennaInfo antenna;
      antenna.timeStart = Sinex::str2time(line, 16, FALSE);
      antenna.timeEnd   = Sinex::str2time(line, 29, TRUE);
      antenna.name      = String::trim(line.substr(42, 15));
      antenna.radome    = String::trim(line.substr(58, 4));
      antenna.serial    = String::trim(line.substr(63, 5));
      if(antenna.radome == "none" || antenna.radome == "NONE")
        antenna.radome = "";
      if(antenna.serial == "n/a" || antenna.serial == "N/A")
        antenna.serial = "";

      station.antenna.push_back(antenna);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

// +SITE/ECCENTRICITY
void GnssStationLog2StationInfo::readSinexEccentricity(std::string &line, GnssStationInfo &station) const
{
  try
  {
    line.resize(std::max(line.size(), UInt(72)), ' ');
    std::string site = String::trim(line.substr(1, 4));
    std::transform(site.begin(), site.end(), site.begin(), ::toupper);

    if(site == station.markerName)
    {
      Time timeStart     = Sinex::str2time(line, 16, FALSE);
      Time timeEnd       = Sinex::str2time(line, 29, TRUE);
      std::string refSys = String::trim(line.substr(42, 3));
      Vector3d eccentricity(String::toDouble(line.substr(55, 8)), String::toDouble(line.substr(64, 8)), String::toDouble(line.substr(46, 8)));

      if(refSys  != "UNE")
      {
        logWarning<<station.markerName<<"."<<station.markerNumber<<": unknown eccentricity reference system '"<<refSys<<"' in sinex file for time period "
                  <<timeStart.dateTimeStr()<<" to "<<timeEnd.dateTimeStr() <<Log::endl;
        return;
      }

      UInt iStart = station.findAntenna(timeStart);
      UInt iEnd   = station.findAntenna(timeEnd - seconds2time(0.1));

      if(iStart == NULLINDEX || iEnd == NULLINDEX)
      {
        logWarning<<station.markerName<<"."<<station.markerNumber<<": no antenna found in sinex file to match eccentricity for time period "
                  <<timeStart.dateTimeStr()<<" to "<<timeEnd.dateTimeStr() <<Log::endl;
        return;
      }

      // split last antenna that fits eccentricity time period if eccentricity ends before antenna
      if(iEnd >= iStart && station.antenna.at(iEnd).timeEnd > timeEnd)
      {
        station.antenna.insert(station.antenna.begin()+iEnd, station.antenna.at(iEnd));
        station.antenna.at(iEnd).timeEnd     = timeEnd;
        station.antenna.at(iEnd+1).timeStart = timeEnd;
      }

      for(UInt i = iStart; i <= iEnd; i++)
        station.antenna.at(i).position = eccentricity;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
