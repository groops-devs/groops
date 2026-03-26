/***********************************************/
/**
* @file gnssTroposphere2TropoSinex.cpp
*
* @brief Convert GNSS troposphere data from GROOPS format to IGS SINEX_TRO format.
*
* @author Barbara Suesser-Rechberger
* @author Sebastian Strasser
* @date 2019-05-29
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Convert GNSS troposphere data from GROOPS format as estimated by \program{GnssProcessing}
to \href{https://files.igs.org/pub/data/format/sinex_tro_v2.00.pdf}{IGS SINEX TRO} format.

For each station folling files are needed:
\begin{itemize}
  \item \configFile{inputfileStationInfo}{platform},
  \item \configFile{inputfileTroposphereData}{instrument}(MISVALUES),
  \item optional standard deviations with \configFile{inputfileTroposphereSigmas}{instrument}(MISVALUES),
  \item \configFile{inputfilePosition}{instrument}(VECTOR3D).
\end{itemize}
The \configFile{inputfileAntennaDefinition}{gnssAntennaDefinition} contains antenna center offsets and variations
of all used antennas. Created via \program{GnssAntex2AntennaDefinition} or \program{GnssAntennaDefinitionCreate}.

For considering the geoid height use \configFile{inputfileGriddedGeoidHeight}{griddedData}
as it might be computed by \program{Gravityfield2GriddedData}.
The height closest to the station's position is used in each case.

See also \program{GnssProcessing}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "files/fileMatrix.h"
#include "files/fileGriddedData.h"
#include "files/fileInstrument.h"
#include "inputOutput/system.h"
#include "inputOutput/file.h"
#include "inputOutput/fileSinex.h"
#include "files/fileStringTable.h"
#include "files/filePlatform.h"

/***** CLASS ***********************************/

/** @brief Convert GNSS troposphere data from GROOPS format to IGS SINEX_TRO format.
* @ingroup programsConversionGroup */
class GnssTroposphere2TropoSinex
{
  static std::string resize(std::string str, UInt length) {str.resize(length, ' '); return str;}

public:
  void run(Config &config, Parallel::CommunicatorPtr comm);

  class Station
  {
  public:
    FileName      fileNameTropoData, fileNameTropoSigmas, fileNameStationInfo, fileNamePosition;
    Time          timesMid;
    Platform      platform;
    Vector3d      position;
    Double        geoidHeight;
    MiscValuesArc tropData, tropSigmas;
  };
};

GROOPS_REGISTER_PROGRAM(GnssTroposphere2TropoSinex, SINGLEPROCESS, "Convert GNSS troposphere data from GROOPS format to SINEX_TRO format.", Conversion, Gnss)

/***********************************************/

template<> Bool readConfig(Config &config, const std::string &name, GnssTroposphere2TropoSinex::Station &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  if(!readConfigSequence(config, name, mustSet, defaultValue, annotation))
    return FALSE;
  readConfig(config, "inputfileStationInfo",       var.fileNameStationInfo, Config::MUSTSET,  "{groopsDataDir}/gnss/receiverStation/stationInfo/stationInfo.{station}.xml", "platform file");
  readConfig(config, "inputfileTroposphereData",   var.fileNameTropoData,   Config::MUSTSET,  "", "Troposphere data estimates (columns: mjd, trodry, trowet, tgndry, tgnwet, tgedry, tgewet)");
  readConfig(config, "inputfileTroposphereSigmas", var.fileNameTropoSigmas, Config::OPTIONAL, "", "Troposphere data sigmas (columns: mjd, sigma_trowet, sigma_tgnwet, sigma_tgewet)");
  readConfig(config, "inputfilePosition",          var.fileNamePosition,    Config::OPTIONAL, "", "Precise station position (columns: mjd, x, y, z [m in TRF])");
  endSequence(config);
  return TRUE;
}

/***********************************************/

void GnssTroposphere2TropoSinex::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName             fileNameTropoSinex;
    std::vector<Station> stations;
    FileName             fileNameAntennaDef, fileNameGeoidHeight;
    std::string          agencyCode, observationCode, solutionContents, gnssSystems, timeSystem, tropoModelingMethod, oceanTideModel, atmosphericTideModel, geoidModel, metDataSource, observationWeighting;
    std::string          description, output, contact, software, hardware, input, versionNumber, systemCode, remark, antennaModel, aPrioriTropoModel, tropoMappingFunction, gradientMappingFunction;
    Time                 timeStart, timeEnd;
    std::vector<std::string> comment;
    FileName             fileNameComment;
    Double               troSampling=NAN_EXPR, dataSampling=NAN_EXPR, elevationCutoff=NAN_EXPR;

    readConfig(config, "outputfileTropoSinex",        fileNameTropoSinex,      Config::MUSTSET,  "",    "");
    readConfig(config, "station",                     stations,                Config::MUSTSET,  "",    "");
    readConfig(config, "inputfileAntennaDefinition",  fileNameAntennaDef,      Config::MUSTSET,  "{groopsDataDir}/gnss/receiverStation/antennaDefinition/igs/igs20/antennaDefinition_igs20.dat", "station phase centers and variations");
    readConfig(config, "inputfileGriddedGeoidHeight", fileNameGeoidHeight,     Config::OPTIONAL, "",    "value closest to the station's position is used in each case");
    readConfig(config, "dataSamplingInterval",        dataSampling,            Config::OPTIONAL, "",    "[sec] GNSS data sampling rate");
    readConfig(config, "tropoSamplingInterval",       troSampling,             Config::OPTIONAL, "",    "[sec] Tropospheric parameter sampling interval");
    readConfig(config, "tropoModelingMethod",         tropoModelingMethod,     Config::OPTIONAL, "Least Squares", "Tropospheric estimation method: Filter, Smoother, Least Squares, Piece-Wise Linear Interpolation");
    readConfig(config, "aPrioriTropoModel",           aPrioriTropoModel,       Config::OPTIONAL, "VMF3 Operational 1x1° gridded + gradients (V3GR)", "A priori tropospheric model used");
    readConfig(config, "tropoMappingFunction",        tropoMappingFunction,    Config::OPTIONAL, "Vienna Mapping Functions 3 (VMF3)", "Name of mapping function used for hydrostatic and wet delay");
    readConfig(config, "gradientMappingFunction",     gradientMappingFunction, Config::OPTIONAL, "Chen & Herring (1997), C_dry=0.0031, C_wet=0.0007", "Name of mapping function used for gradients");
    readConfig(config, "metDataSource",               metDataSource,           Config::OPTIONAL, "",    "source of surface meteorological observations used (see format desc.)");
    readConfig(config, "observationWeighting",        observationWeighting,    Config::OPTIONAL, "",    "observation weighting model applied");
    readConfig(config, "elevationCutoff",             elevationCutoff,         Config::OPTIONAL, "5",   "[deg]");
    readConfig(config, "gnssSystems",                 gnssSystems,             Config::OPTIONAL, "",    "G=GPS, R=GLONASS, E=Galileo, C=BeiDou");
    readConfig(config, "timeSystem",                  timeSystem,              Config::OPTIONAL, "G",   "G (GPS) or UTC");
    readConfig(config, "oceanTideModel",              oceanTideModel,          Config::OPTIONAL, "",    "Name of applied Ocean tide loading model");
    readConfig(config, "atmosphericTideModel",        atmosphericTideModel,    Config::OPTIONAL, "",    "Name of applied Atmospheric tide loading model");
    readConfig(config, "geoidModel",                  geoidModel,              Config::OPTIONAL, "",    "Geoid model name for undulation values");
    readConfig(config, "systemCode",                  systemCode,              Config::OPTIONAL, "",    "Terrestrial reference system code");
    readConfig(config, "remark",                      remark,                  Config::OPTIONAL, "TUG", "Remark used to identify the origin of the coordinates (AC acronym)");
    readConfig(config, "antennaCalibrationModel",     antennaModel,            Config::OPTIONAL, "",    "e.g. IGS20_WWWW (WWWW = ANTEX release GPS week)");
    if(readConfigSequence(config, "sinexTroHeader", Config::MUSTSET, "", ""))
    {
      readConfig(config, "agencyCode",       agencyCode,       Config::OPTIONAL, "TUG", "Identify the agency providing the data");
      readConfig(config, "timeStart",        timeStart,        Config::OPTIONAL, "",    "Start time of the data");
      readConfig(config, "timeEnd",          timeEnd,          Config::OPTIONAL, "",    "End time of the data ");
      readConfig(config, "observationCode",  observationCode,  Config::OPTIONAL, "P",   "Technique used to generate the SINEX solution");
      readConfig(config, "solutionContents", solutionContents, Config::OPTIONAL, "MIX", "Marker name if single station, MIX if multiple stations");
      readConfig(config, "description",      description,      Config::OPTIONAL, "",    "Organizitions gathering/alerting the file contents");
      readConfig(config, "output",           output,           Config::OPTIONAL, "",    "Description of the file contents");
      readConfig(config, "contact",          contact,          Config::OPTIONAL, "",    "Address of the relevant contact e-mail");
      readConfig(config, "software",         software,         Config::OPTIONAL, "GROOPS (https://github.com/groops-devs/groops)", "Software used to generate the file");
      readConfig(config, "hardware",         hardware,         Config::OPTIONAL, "",    "Computer hardware on which above software was run");
      readConfig(config, "input",            input,            Config::OPTIONAL, "",    "Brief description of the input used to generate this solution");
      readConfig(config, "versionNumber",    versionNumber,    Config::OPTIONAL, "",    "Unique identifier of the product, same as in file name, e.g. 000");
      readConfig(config, "inputfileComment", fileNameComment,  Config::OPTIONAL, "",    "comments in the comment block from a file (truncated at 80 characters per line)");
      readConfig(config, "comment",          comment,          Config::OPTIONAL, "",    "comments in the comment block");
      endSequence(config);
    }
    if(isCreateSchema(config)) return;

    // ==================================================

    logStatus<<"reading station antenna definitions from <"<<fileNameAntennaDef<<">"<<Log::endl;
    std::vector<GnssAntennaDefinitionPtr> antennaDefinitionList;
    readFileGnssAntennaDefinition(fileNameAntennaDef, antennaDefinitionList);

    GriddedData gridGeoid;
    if(!fileNameGeoidHeight.empty())
    {
      logStatus<<"reading geoid heights <"<<fileNameGeoidHeight<<">"<<Log::endl;
      readFileGriddedData(fileNameGeoidHeight, gridGeoid);
    }

    logStatus<<"reading station infos"<<Log::endl;
    for(auto &station : stations)
    {
      station.tropData   = InstrumentFile::read(station.fileNameTropoData);
      station.tropSigmas = InstrumentFile::read(station.fileNameTropoSigmas);
      Arc::checkSynchronized({station.tropData, station.tropSigmas});
      station.timesMid = 0.5*(station.tropData.front().time+station.tropData.back().time);

      readFilePlatform(station.fileNameStationInfo, station.platform);
      station.platform.fillGnssAntennaDefinition(antennaDefinitionList);

      // position
      station.position = station.platform.approxPosition;
      if(!station.fileNamePosition.empty())
      {
        Vector3dArc arc = InstrumentFile::read(station.fileNamePosition);
        auto iter = std::min_element(arc.begin(), arc.end(), [&](const Epoch &e1, const Epoch &e2)
                                    {return std::fabs((e1.time-station.timesMid).mjd()) < std::fabs((e2.time-station.timesMid).mjd());});
        if(!arc.size() || ((arc.size() > 1) && (std::fabs((iter->time-station.timesMid).mjd()) > 0.5*medianSampling(arc.times()).mjd())))
          throw(Exception("No position found in <"+station.fileNamePosition.str()+">"));
        station.position = iter->vector3d;
      }

      // find nearest geoid height
      station.geoidHeight = NAN_EXPR;
      if(gridGeoid.points.size() && gridGeoid.values.size())
      {
        auto iter = std::min_element(gridGeoid.points.begin(), gridGeoid.points.end(), [&](const Vector3d &p1, const Vector3d &p2)
                                     {return (p1-station.position).quadsum() < (p2-station.position).quadsum();});
        station.geoidHeight = gridGeoid.values.front().at(std::distance(gridGeoid.points.begin(), iter));
      }
    }

    // ==================================================

    // write SINEX_TRO file
    // ----------------
    logStatus<<"write SINEX_TRO file <"<<fileNameTropoSinex<<">"<<Log::endl;
    Sinex sinex;
    // SINEX_TRO header line
    std::stringstream ss;
    ss<<"%=TRO 2.00 "<<std::setw(3)<<agencyCode.substr(0, 3)<<" "<<Sinex::time2str(System::now(), TRUE)<<" "<<std::setw(3)<<agencyCode.substr(0, 3);
    ss<<" "<<Sinex::time2str(timeStart, TRUE)<<" "<<Sinex::time2str(timeEnd, TRUE)<<" "<<std::setw(1)<<observationCode.substr(0, 1);
    ss<<" "<<solutionContents.substr(0, 4);
    sinex.header = ss.str();
    sinex.footer = "%=ENDTRO";

    std::vector<std::string> columns = {"TROTOT", "TGNTOT", "TGETOT"};
    if(std::any_of(stations.begin(), stations.end(), [](auto &station) {return station.tropSigmas.size();}))
      columns = {"TROTOT", "STDDEV", "TGNTOT", "STDDEV", "TGETOT", "STDDEV"};

    // Block: FILE/REFERENCE
    {
      SinexBlockPtr block = sinex.addBlock("FILE/REFERENCE");
      *block<<"*INFO_TYPE_________ INFO________________________________________________________"<<std::endl;
      if(!description.empty())   *block<<" DESCRIPTION        "<<description<<std::endl;
      if(!output.empty())        *block<<" OUTPUT             "<<output<<std::endl;
      if(!contact.empty())       *block<<" CONTACT            "<<contact<<std::endl;
      if(!software.empty())      *block<<" SOFTWARE           "<<software<<std::endl;
      if(!hardware.empty())      *block<<" HARDWARE           "<<hardware<<std::endl;
      if(!input.empty())         *block<<" INPUT              "<<input<<std::endl;
      if(!versionNumber.empty()) *block<<" VERSION NUMBER     "<<versionNumber<<std::endl;
    }

    // Block: TROP/DESCRIPTION
    {
      SinexBlockPtr block = sinex.addBlock("TROP/DESCRIPTION");
      *block<<"*_________KEYWORD____________ ___VALUES(S)______________________________________"<<std::endl;
      *block<<" TROPO PARAMETER NAMES       "; for(const auto &c : columns)         {*block<<" "<<c;}    *block<<std::endl;
      *block<<" TROPO PARAMETER UNITS       "; for(UInt i=0; i<columns.size(); i++) {*block<<"  1e+03";} *block<<std::endl;
      *block<<" TROPO PARAMETER WIDTH       "; for(UInt i=0; i<columns.size(); i++) {*block<<"      6";} *block<<std::endl;
      if(!tropoModelingMethod.empty())     *block<<" TROPO MODELING METHOD        "<<tropoModelingMethod<<std::endl;
      if(!std::isnan(troSampling))         *block<<" TROPO SAMPLING INTERVAL      "<<troSampling%"%f"s<<std::endl;
      if(!metDataSource.empty())           *block<<" SOURCE OF MET/DATA           "<<metDataSource<<std::endl;
      if(!aPrioriTropoModel.empty())       *block<<" A PRIORI TROPOSPHERE         "<<aPrioriTropoModel<<std::endl;
      if(!tropoMappingFunction.empty())    *block<<" TROPO MAPPING FUNCTION       "<<tropoMappingFunction<<std::endl;
      if(!gradientMappingFunction.empty()) *block<<" GRADS MAPPING FUNCTION       "<<gradientMappingFunction<<std::endl;
      if(!std::isnan(dataSampling))        *block<<" DATA SAMPLING INTERVAL       "<<dataSampling%"%f"s<<std::endl;
      if(!std::isnan(elevationCutoff))     *block<<" ELEVATION CUTOFF ANGLE       "<<elevationCutoff%"%f"s<<std::endl;
      if(!observationWeighting.empty())    *block<<" OBSERVATION WEIGHTING        "<<observationWeighting<<std::endl;
      if(!gnssSystems.empty())             *block<<" GNSS SYSTEMS                 "<<gnssSystems<<std::endl;
      if(!aPrioriTropoModel.empty())       *block<<" TIME SYSTEM                  "<<timeSystem<<std::endl;
      if(!timeSystem.empty())              *block<<" OCEAN TIDE LOADING MODEL     "<<oceanTideModel<<std::endl;
      if(!atmosphericTideModel.empty())    *block<<" ATMOSPH TIDE LOADING MODEL   "<<atmosphericTideModel<<std::endl;
      if(!geoidModel.empty())              *block<<" GEOID MODEL                  "<<geoidModel<<std::endl;
    }

    // Block: FILE/COMMENT
    if(comment.size() || !fileNameComment.empty())
    {
      SinexBlockPtr block = sinex.addBlock("SITE/ID");
      if(!fileNameComment.empty())
      {
        InFile commentFile(fileNameComment);
        std::string line;
        while(std::getline(commentFile, line))
          *block<<" "<<line<<std::endl;
      }
      for(const auto &line : comment)
        *block<<" "<<line<<std::endl;
    }

    // Block: SITE/ID
    {
      SinexBlockPtr block = sinex.addBlock("SITE/ID");
      *block<<"*STATION__ PT __DOMES__ T _STATION_DESCRIPTION__ _LONGITUDE _LATITUDE_ _HGT_ELI_ _HGT_MSL_"<<std::endl;
      Ellipsoid ellipsoid;
      for(const auto &station : stations)
      {
        Angle lon, lat;
        Double ellipsoidHeight;
        ellipsoid(station.position, lon, lat, ellipsoidHeight);
        *block<<" "<<String::upperCase(resize(station.platform.markerName, 9))<<"  A "<<resize(station.platform.markerNumber, 9)<<" P "<<resize(station.platform.comment, 22)<<" "
              <<(std::fmod(lon+2*PI, 2*PI)*RAD2DEG)%"%10.6f "s<<(lat*RAD2DEG)%"%10.6f "s<<ellipsoidHeight%"%9.3f "s<<(std::isnan(station.geoidHeight) ? "" : station.geoidHeight % "%9.3f"s)<<std::endl;
      }
    }

    // Block: SITE/COORDINATES
    {
      SinexBlockPtr block = sinex.addBlock("SITE/COORDINATES");
      *block<<"*STATION__ PT SOLN T __DATA_START__ __DATA_END____ __STA_X_____  __STA_Y_____  __STA_Z_____ SYSTEM REMRK"<<std::endl;
      for(const auto &station : stations)
      {
        *block<<" "<<String::upperCase(resize(station.platform.markerName, 9))<<"  A    1 P "<<Sinex::time2str(timeStart, TRUE)<<" "<<Sinex::time2str(timeEnd, TRUE)<<" "<<station.position.x() % "%12.3f "s;
        *block<<" "<<station.position.y()%"%12.3f "s<<" "<<station.position.z()%"%12.3f "s<<resize(systemCode, 6)<<" "<<resize(remark, 5)<<std::endl;
      }
    }

    // BLOCK: SITE/RECEIVER
    {
      SinexBlockPtr block = sinex.addBlock("SITE/RECEIVER");
      *block<<"*STATION__ PT SOLN T __DATA_START__ __DATA_END____ DESCRIPTION_________ S/N_________________ FIRMW______"<<std::endl;
      for(const auto &station : stations)
      {
        auto recv = station.platform.findEquipment<PlatformGnssReceiver>(station.timesMid);
        if(!recv)
        {
          logWarning<<station.platform.markerName<<": no receiver found at "<<station.timesMid.dateTimeStr()<<Log::endl;
          continue;
        }
        *block<<" "<<String::upperCase(resize(station.platform.markerName, 9))<<"  A    1 P "<<Sinex::time2str(recv->timeStart, TRUE)<<" "<<Sinex::time2str(recv->timeEnd, TRUE)<<" "<<resize(recv->name, 20)
              <<" "<<resize(recv->serial.empty() ? std::string(20, '-') : recv->serial, 20)<<" "<<resize(recv->version.empty() ? std::string(11, '-') : recv->version, 11)<<std::endl;
      }
    }

    // BLOCK SITE/ANTENNA
    {
      SinexBlockPtr block = sinex.addBlock("SITE/ANTENNA");
      *block<<"*STATION__ PT SOLN T __DATA_START__ __DATA_END____ DESCRIPTION_________ S/N_________________ PCV_MODEL_"<<std::endl;
      for(const auto &station : stations)
      {
        auto ant = station.platform.findEquipment<PlatformGnssAntenna>(station.timesMid);
        if(!ant)
        {
          logWarning<<station.platform.markerName<<": no antenna found at "<<station.timesMid.dateTimeStr()<<Log::endl;
          continue;
        }
        *block<<" "<<String::upperCase(resize(station.platform.markerName, 9))<<"  A    1 P "<<Sinex::time2str(ant->timeStart, TRUE)<<" "<<Sinex::time2str(ant->timeEnd, TRUE)<<" "
              <<resize(ant->name, 15)<<" "<<(ant->radome.empty() ? "NONE" : resize(ant->radome, 4))<<" "<<resize(ant->serial.empty() ? std::string(20, '-') : ant->serial, 20)<<" "
              <<resize(antennaModel, 10)<<std::endl;
      }
    }

    // BLOCK SITE/ECCENTRICITY
    {
      SinexBlockPtr block = sinex.addBlock("SITE/ECCENTRICITY");
      *block<<"*"<<std::setw(80)<<"UP______ NORTH___ EAST____"<<std::endl;
      *block<<"*STATION__ PT SOLN T __DATA_START__ __DATA_END____ AXE ARP->BENCHMARK(M)_________"<<std::endl;
      for(const auto &station : stations)
      {
        auto ant = station.platform.findEquipment<PlatformGnssAntenna>(station.timesMid);
        if(!ant)
        {
          logWarning<<station.platform.markerName<<": no antenna found at "<<station.timesMid.dateTimeStr()<<Log::endl;
          continue;
        }
        *block<<" "<<String::upperCase(resize(station.platform.markerName, 9))<<"  A    1 P "<<Sinex::time2str(ant->timeStart, TRUE)<<" "
              <<Sinex::time2str(ant->timeEnd, TRUE)<<" UNE "<<ant->position.z() % "%8.4f "s<<ant->position.x() % "%8.4f "s<<ant->position.y() % "%8.4f"s<<std::endl;
      }
    }

    // BLOCK: TROP/SOLUTION
    {
      SinexBlockPtr block = sinex.addBlock("TROP/SOLUTION");
      *block<<"*STATION__ ____EPOCH_____"; for(const auto &c : columns) {*block<<" "<<c;}; *block<<std::endl;
      for(const auto &station : stations)
        for(UInt idEpoch=0; idEpoch<station.tropData.size(); idEpoch++)
        {
          const Vector data = 1e3 * station.tropData.at(idEpoch).values;
          *block<<" "<<String::upperCase(resize(station.platform.markerName, 9))<<" "<<Sinex::time2str(station.tropData.at(idEpoch).time, TRUE);
          if(!station.tropSigmas.size())
            *block<<(data(0)+data(1))%" %6.1f"s<<(data(2)+data(3))%" %6.2f"s<<(data(4)+data(5))%" %6.2f"s<<std::endl;
          else
          {
            const Vector sigmas = 1e3 * station.tropSigmas.at(idEpoch).values;
            *block<<(data(0)+data(1))%" %6.1f"s<<sigmas(0)%" %6.1f"s
                  <<(data(2)+data(3))%" %6.2f"s<<sigmas(1)%" %6.2f"s
                  <<(data(4)+data(5))%" %6.2f"s<<sigmas(2)%" %6.2f"s<<std::endl;
          }
        }
    }

    writeFileSinex(fileNameTropoSinex, sinex);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
