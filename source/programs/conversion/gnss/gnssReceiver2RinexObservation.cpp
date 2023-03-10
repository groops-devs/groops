/***********************************************/
/**
* @file gnssReceiver2RinexObservation.cpp
*
* @brief Converts GnssReceiver Instrument file to RINEX.
*
* @author Torsten Mayer-Guerr
* @date 2023-03-05
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Converts a \configFile{inputfileGnssReceiver}{instrument} into a
\href{https://files.igs.org/pub/data/format/rinex_4.00.pdf}{RINEX} observation file.
The \configFile{inputfileStationInfo}{platform} contains the antenna
and receiver information for the RINEX header.
The \configClass{useType}{gnssType} and \configClass{ignoreType}{gnssType} can be used to filter
the observation types that will be exported.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "inputOutput/file.h"
#include "inputOutput/system.h"
#include "files/filePlatform.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Converts GnssReceiver Instrument file to RINEX.
* @ingroup programsConversionGroup */
class GnssReceiver2RinexObservation
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GnssReceiver2RinexObservation, SINGLEPROCESS, "Converts GnssReceiver Instrument file to RINEX.", Conversion, Gnss, Instrument)

/***********************************************/

void GnssReceiver2RinexObservation::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameOut, fileNameObs, fileNamePlatform;
    std::string observer, angency;
    std::vector<std::string> comments;
    std::vector<GnssType> useType, ignoreType;

    readConfig(config, "outputfileRinexObservation", fileNameOut,      Config::MUSTSET,  "",    "RINEX observation file");
    readConfig(config, "inputfileGnssReceiver",      fileNameObs,      Config::MUSTSET,  "",    "GNSS instrument file");
    readConfig(config, "inputfileStationInfo",       fileNamePlatform, Config::MUSTSET,  "",    "antenna and receiver info");
    readConfig(config, "comment",                    comments,         Config::OPTIONAL, "",    "write comments at begin of header");
    readConfig(config, "observer",                   observer,         Config::OPTIONAL, "TUG", "header information");
    readConfig(config, "angency",                    angency,          Config::OPTIONAL, "TUG", "header information");
    readConfig(config, "useType",                    useType,          Config::OPTIONAL, "",    "only use observations that match any of these patterns");
    readConfig(config, "ignoreType",                 ignoreType,       Config::OPTIONAL, "",    "ignore observations that match any of these patterns");
    if(isCreateSchema(config)) return;

    logStatus<<"read instrument file <"<<fileNameObs<<">"<<Log::endl;
    GnssReceiverArc arc = InstrumentFile::read(fileNameObs);
    if(!arc.size())
    {
      logWarning<<"Empty observation file"<<Log::endl;
      return;
    }

    // find all types
    std::map<GnssType, std::set<GnssType>> types; // types per system
    std::set<GnssType> glonass;                   // GLONASS frequencyNumbers
    for(auto &epoch : arc)
      for(GnssType satType : epoch.satellite)
        for(UInt idType = GnssType::index(epoch.obsType, satType); (idType < epoch.obsType.size()) && (epoch.obsType.at(idType) == satType); idType++)
          if((!useType.size() || (epoch.obsType.at(idType)+satType).isInList(useType)) && !(epoch.obsType.at(idType)+satType).isInList(ignoreType))
          {
            types[satType & GnssType::SYSTEM].insert(epoch.obsType.at(idType));
            if(satType == GnssType::GLONASS)
              glonass.insert(satType);
          }

    logStatus<<"read platform file <"<<fileNamePlatform<<">"<<Log::endl;
    Platform platform;
    readFilePlatform(fileNamePlatform, platform);

    const Time timeStart = arc.front().time;
    Vector3d com  = platform.referencePoint(timeStart);
    auto antenna  = platform.findEquipment<PlatformGnssAntenna>(timeStart);
    if(!antenna)
      throw(Exception(platform.markerName+"."+platform.markerNumber+": no antenna definition found at "+timeStart.dateTimeStr()));
    auto receiver = platform.findEquipment<PlatformGnssReceiver>(timeStart);

    auto fixedSize = [](const std::string &str, UInt size)
    {
      return (str+std::string(size, ' ')).substr(0, size);
    };

    // header
    // ------
    logStatus<<"write RINEX file <"<<fileNameOut<<">"<<Log::endl;
    OutFile file(fileNameOut);
    file<<"     4.00           O                   M                   RINEX VERSION / TYPE"<<std::endl;
    for(auto &comment : comments)
      file<<fixedSize(comment, 60)                                  <<"COMMENT             "<<std::endl;
    file<<"GROOPS              "<<fixedSize(angency, 20)<<System::now()%"%y%m%d %H%M%S LCL PGM / RUN BY / DATE"s<<std::endl;
    file<<fixedSize(observer, 20)<<fixedSize(angency, 40)           <<"OBSERVER / AGENCY   "<<std::endl;
    file<<fixedSize(platform.markerName, 60)                        <<"MARKER NAME         "<<std::endl;
    if(!platform.markerNumber.empty())
      file<<fixedSize(platform.markerNumber, 20)<<"                                        MARKER NUMBER       "<<std::endl;
    if(receiver)
      file<<fixedSize(receiver->serial, 20)<<fixedSize(receiver->name, 20)<<fixedSize(receiver->version, 20)<<"REC # / TYPE / VERS "<<std::endl;
    file<<fixedSize(antenna->serial, 20)<<fixedSize(antenna->name, 16)<<fixedSize(antenna->radome, 4)<<"                    ANT # / TYPE        "<<std::endl;
    file<<platform.approxPosition.x()%"%14.4f"s<<platform.approxPosition.y()%"%14.4f"s<<platform.approxPosition.z()%"%14.4f                  APPROX POSITION XYZ "s<<std::endl;
    file<<antenna->position.z()%"%14.4f"s<<antenna->position.y()%"%14.4f"s<<antenna->position.x()%"%14.4f                  ANTENNA: DELTA H/E/N"s<<std::endl;
    file<<antenna->position.x()%"%14.4f"s<<antenna->position.y()%"%14.4f"s<<antenna->position.z()%"%14.4f                  ANTENNA: DELTA X/Y/Z"s<<std::endl;
    const Vector3d boresight = antenna->local2antennaFrame.inverseTransform(Vector3d(0,0,1));
    const Vector3d zerodir   = antenna->local2antennaFrame.inverseTransform(Vector3d(1,0,0));
    file<<boresight.x()%"%14.4f"s<<boresight.y()%"%14.4f"s<<boresight.z()%"%14.4f                  ANTENNA: B.SIGHT XYZ"s<<std::endl;
    file<<zerodir.x()  %"%14.4f"s<<zerodir.y()  %"%14.4f"s<<zerodir.z()  %"%14.4f                  ANTENNA: ZERODIR XYZ"s<<std::endl;
    if(com.r() >= 1e-4)
      file<<com.x()%"%14.4f"s<<com.y()%"%14.4f"s<<com.z()%"%14.4f                  CENTER OF MASS: XYZ "s<<std::endl;
    for(auto &system : types)
    {
      file<<system.first.str().at(3)<<system.second.size()%"  %3i"s;
      UInt count = 0;
      for(auto &type : system.second)
      {
        if(++count > 13) // continuation line?
        {
          count = 1;
          file<<"  SYS / # / OBS TYPES "<<std::endl<<"      ";
        }
        file<<" "<<type.str().substr(0, 3);
      }
      file<<std::string(60-6-4*count, ' ')<<"SYS / # / OBS TYPES "<<std::endl;
    }
    file<<timeStart      %"  %y    %m    %d    %H    %M  %11.7S     GPS         TIME OF FIRST OBS "s<<std::endl;
    file<<arc.back().time%"  %y    %m    %d    %H    %M  %11.7S     GPS         TIME OF LAST OBS  "s<<std::endl;
    if(glonass.size())
    {
      file<<glonass.size()%"%3i "s;
      UInt count = 0;
      for(auto &type : glonass)
      {
        if(++count > 8) // continuation line?
        {
          count = 1;
          file<<"GLONASS SLOT / FRQ #"<<std::endl<<"    ";
        }
        file<<type.str().substr(3, 3)<<type.frequencyNumber()%" %2i "s;
      }
      file<<std::string(60-4-7*count, ' ')<<"GLONASS SLOT / FRQ #"<<std::endl;
    }
    file<<"                                                            END OF HEADER       "<<std::endl;

    // data section
    // ------------
    for(auto &epoch : arc)
    {
      // antenna change?
      if(antenna != platform.findEquipment<PlatformGnssAntenna>(epoch.time))
      {
        antenna = platform.findEquipment<PlatformGnssAntenna>(epoch.time);
        if(!antenna)
          throw(Exception(platform.markerName+"."+platform.markerNumber+": no antenna definition found at "+epoch.time.dateTimeStr()));
        file<<"> "<<antenna->timeStart%"%y %m %d %H %M%11.7S  4  5"s<<std::endl;
        file<<fixedSize(antenna->serial, 20)<<fixedSize(antenna->name, 16)<<fixedSize(antenna->radome, 4)<<"                    ANT # / TYPE        "<<std::endl;
        file<<antenna->position.z()%"%14.4f"s<<antenna->position.y()%"%14.4f"s<<antenna->position.x()%"%14.4f                  ANTENNA: DELTA H/E/N"s<<std::endl;
        file<<antenna->position.x()%"%14.4f"s<<antenna->position.y()%"%14.4f"s<<antenna->position.z()%"%14.4f                  ANTENNA: DELTA X/Y/Z"s<<std::endl;
        const Vector3d boresight = antenna->local2antennaFrame.inverseTransform(Vector3d(0,0,1));
        const Vector3d zerodir   = antenna->local2antennaFrame.inverseTransform(Vector3d(1,0,0));
        file<<boresight.x()%"%14.4f"s<<boresight.y()%"%14.4f"s<<boresight.z()%"%14.4f                  ANTENNA: B.SIGHT XYZ"s<<std::endl;
        file<<zerodir.x()  %"%14.4f"s<<zerodir.y()  %"%14.4f"s<<zerodir.z()  %"%14.4f                  ANTENNA: ZERODIR XYZ"s<<std::endl;
      }

      // receiver change?
      if(receiver != platform.findEquipment<PlatformGnssReceiver>(epoch.time))
      {
        receiver = platform.findEquipment<PlatformGnssReceiver>(epoch.time);
        if(receiver)
        {
          file<<"> "<<receiver->timeStart%"%y %m %d %H %M%11.7S  4  1"s<<std::endl;
          file<<fixedSize(receiver->serial, 20)<<fixedSize(receiver->name, 20)<<fixedSize(receiver->version, 20)<<"REC # / TYPE / VERS "<<std::endl;
        }
      }

      // center of mass change?
      if((com-platform.referencePoint(epoch.time)).r() >= 1e-4)
      {
        com = platform.referencePoint(epoch.time);
        file<<"> "<<epoch.time%"%y %m %d %H %M%11.7S  4  1"s<<std::endl;
        file<<com.x()%"%14.4f"s<<com.y()%"%14.4f"s<<com.z()%"%14.4f                  CENTER OF MASS: XYZ "s<<std::endl;
      }

      // Has satellite valid observation types?
      auto isSat = [&](GnssType satType)
      {
        for(UInt idType = GnssType::index(epoch.obsType, satType); (idType < epoch.obsType.size()) && (epoch.obsType.at(idType) == satType); idType++)
          if((!useType.size() || (epoch.obsType.at(idType)+satType).isInList(useType)) && !(epoch.obsType.at(idType)+satType).isInList(ignoreType))
            return TRUE;
        return FALSE;
      };

      const UInt countSat = std::count_if(epoch.satellite.begin(), epoch.satellite.end(), isSat);
      if(!countSat)
        continue;
      file<<"> "<<epoch.time%"%y %m %d %H %M%11.7S  0"s<<countSat%"%3i"s;
      if(epoch.clockError)
        file<<epoch.clockError%"      %15.12f"s;
      file<<std::endl;

      UInt idObs = 0, idx;
      for(GnssType satType : epoch.satellite)
      {
        UInt idType = GnssType::index(epoch.obsType, satType); // find type for the satellite system
        if(isSat(satType))
        {
          file<<satType.str().substr(3, 3);
          for(GnssType type : types.at(satType & GnssType::SYSTEM))
            if(type.isInList(epoch.obsType, idx) && (!useType.size() || (type+satType).isInList(useType)) && !(type+satType).isInList(ignoreType))
              file<<(epoch.observation.at(idObs+idx-idType) / ((type == GnssType::PHASE) ? (type+satType).wavelength() : 1.))%"%14.3f  "s; // convert to cycles
            else
              file<<"                ";
          file<<std::endl;
        }

        // to next satellite
        while((idType < epoch.obsType.size()) && (epoch.obsType.at(idType) == satType))
          idType++, idObs++;
      } // for(satType)
    } // for(epoch)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
