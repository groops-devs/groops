/***********************************************/
/**
* @file gnssTroposphere2TropoSinex.cpp
*
* @brief Convert GNSS troposphere data from GROOPS format to IGS SINEX\_TRO format.
*
* @author Barbara Suesser-Rechberger
* @author Sebastian Strasser
* @date 2019-05-29
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Convert GNSS troposphere data from GROOPS format to \href{http://twg.igs.org/documents/sinex_tro_v2.00.pdf}{IGS SINEX TRO} format.

Specification of the station list is done via \config{inputfileStationList}.
\config{inputfileTroposphereData} needs the troposphere data provided from \program{GnssProcessing}.
Additional following station metadata are required: \config{inputfileStationInfo} using file type \file{Station info}{gnssStationInfo},
\config{inputfileAntennaDefinition} using file type \file{Antenna definition}{gnssAntennaDefinition} and \config{inputfileGridPos} which
uses the stations positions provided from \program{GnssProcessing}.
For considering the geoid height use \config{inputfileGeoidHeight}. The geoid height is provided by \program{Gravityfield2GriddedData}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "files/fileMatrix.h"
#include "files/fileGriddedData.h"
#include "inputOutput/system.h"
#include "inputOutput/file.h"
#include "inputOutput/fileSinex.h"
#include "files/fileStringTable.h"
#include "files/fileGnssStationInfo.h"

/***** CLASS ***********************************/

/** @brief Convert GNSS troposphere data from GROOPS format to IGS SINEX\_TRO format.
* @ingroup programsConversionGroup */
class GnssTroposphere2TropoSinex
{
  static std::string resize(std::string str, UInt length) { str.resize(length, ' '); return str; }
  static std::string time2str(Time time) { return Sinex::time2str(time, TRUE); }
  std::string SiteId(const std::vector<GnssStationInfo> &stationInfos, const GriddedData &griddedData, const std::vector <Double> &geoidHeightVector);
  std::string SiteReceiver(const std::vector<GnssStationInfo> &stationInfos, const Time &timeRef);
  std::string TropSolution(const std::vector<GnssStationInfo> &stationInfos, const std::vector<Matrix> &tropDataMatrix, const std::vector<Matrix> &tropSigmasMatrix);
  std::string SiteCoordinates(const std::vector<GnssStationInfo> &stationInfos, const GriddedData &griddedData, const Time &timeStart,
                              const Time &timeEnd, std::string &systemCode, std::string &remark);
  std::string SiteAntenna(const std::vector<GnssStationInfo> &stationInfos, const Time &timeRef, std::string &antennaModel);
  std::string SiteEccentricity(const std::vector<GnssStationInfo> &stationInfos, const Time &timeRef);

public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(GnssTroposphere2TropoSinex, SINGLEPROCESS, "Convert GNSS troposphere data from GROOPS format to SINEX_TRO format.", Conversion, Gnss)

/***********************************************/

void GnssTroposphere2TropoSinex::run(Config &config)
{
  try
  {
    FileName fileNameTropoSinex, fileNameTropoData, fileNameTropoSigmas, fileNameStationList, fileNameStationInfo, fileNameAntennaDef, variableStationName, fileNameGeoidHeight, fileNameGridPos;
    Time timeRef, timeStartObs, timeEndObs, timeStart, timeEnd;
    std::string agencyCode, observationCode, solutionContents, gnssSystems, timeSystem, tropoModelingMethod, oceanTideModel, atmosphericTideModel, geoidModel, metDataSource, observationWeighting;
    std::string description, output, contact, software, hardware, input, versionNumber, blockText, systemCode, remark, antennaModel, aPrioriTropoModel, tropoMappingFunction, gradientMappingFunction;
    std::vector<std::string> comment;
    FileName fileNameComment;
    Double troSampling, dataSampling, elevationCutoff;

    readConfig(config, "outputfileTropoSinex",     fileNameTropoSinex, Config::MUSTSET, "", "");
    if(readConfigSequence(config, "stations", Config::MUSTSET, "", ""))
    {
      readConfig(config, "inputfileStationList",       fileNameStationList, Config::MUSTSET,  "", "ASCII file with station names");
      readConfig(config, "inputfileTroposphereData",   fileNameTropoData,   Config::MUSTSET,  "", "Troposphere data estimates template (columns: mjd, trodry, trowet, tgndry, tgnwet, tgedry, tgewet)");
      readConfig(config, "inputfileTroposphereSigmas", fileNameTropoSigmas, Config::OPTIONAL, "", "Troposphere data sigmas template (columns: mjd, sigma_trowet, sigma_tgnwet, sigma_tgewet)");
      readConfig(config, "inputfileStationInfo",       fileNameStationInfo, Config::MUSTSET,  "{groopsDataDir}/gnss/receiverStation/stationInfo/igs/stationInfo.{station}.xml", "station info file template");
      readConfig(config, "inputfileGeoidHeight",       fileNameGeoidHeight, Config::OPTIONAL, "", "File including geoid height");
      readConfig(config, "inputfileGridPos",           fileNameGridPos,     Config::MUSTSET,  "", "File including stations positions");
      readConfig(config, "inputfileAntennaDefinition", fileNameAntennaDef,  Config::MUSTSET,  "{groopsDataDir}/gnss/receiverStation/antennaDefinition/igs/igs14/antennaDefinition_igs14.dat", "station phase centers and variations (ANTEX)");
      readConfig(config, "variableStationName",        variableStationName, Config::DEFAULT,  "station", "Loop variable for station names from station list");
      readConfig(config, "observationTimeStart",       timeStartObs,        Config::MUSTSET,  "", "Start time for which solution has observations");
      readConfig(config, "observationTimeEnd",         timeEndObs,          Config::MUSTSET,  "", "End time for which solution has observations");
      endSequence(config);
    }
    readConfig(config, "dataSamplingInterval",    dataSampling,            Config::MUSTSET,  "", "[sec] GNSS data sampling rate");
    readConfig(config, "tropoSamplingInterval",   troSampling,             Config::MUSTSET,  "", "[sec] Tropospheric parameter sampling interval");
    readConfig(config, "tropoModelingMethod",     tropoModelingMethod,     Config::MUSTSET,  "", "Tropospheric estimation method: Filter, Smoother, Least Squares, Piece-Wise Linear Interpolation");
    readConfig(config, "aPrioriTropoModel",       aPrioriTropoModel,       Config::MUSTSET,  "", "A priori tropospheric model used");
    readConfig(config, "tropoMappingFunction",    tropoMappingFunction,    Config::MUSTSET,  "", "Name of mapping function used for hydrostatic and wet delay");
    readConfig(config, "gradientMappingFunction", gradientMappingFunction, Config::MUSTSET,  "", "Name of mapping function used for gradients");
    readConfig(config, "metDataSource",           metDataSource,           Config::OPTIONAL, "", "source of surface meteorological observations used (see format desc.)");
    readConfig(config, "observationWeighting",    observationWeighting,    Config::MUSTSET,  "", "observation weighting model applied");
    readConfig(config, "elevationCutoff",         elevationCutoff,         Config::MUSTSET,  "", "[deg]");
    readConfig(config, "gnssSystems",             gnssSystems,             Config::MUSTSET,  "", "G=GPS, R=GLONASS, E=Galileo, C=BeiDou");
    readConfig(config, "timeSystem",              timeSystem,              Config::MUSTSET,  "G", "G (GPS) or UTC");
    readConfig(config, "oceanTideModel",          oceanTideModel,          Config::OPTIONAL, "", "Ocean tide loading model applied");
    readConfig(config, "atmosphericTideModel",    atmosphericTideModel,    Config::OPTIONAL, "", "Atmospheric tide loading model applied");
    readConfig(config, "geoidModel",              geoidModel,              Config::OPTIONAL, "", "Geoid model name for undulation values");
    readConfig(config, "systemCode",              systemCode,              Config::MUSTSET,  "",    "Terrestrial reference system code");
    readConfig(config, "remark",                  remark,                  Config::MUSTSET,  "TUG", "Remark used to identify the origin of the coordinates (AC acronym)");
    readConfig(config, "antennaCalibrationModel", antennaModel,            Config::MUSTSET,  "",    "e.g. IGS14_WWWW (WWWW = ANTEX release GPS week)");
    if(readConfigSequence(config, "sinexTroHeader", Config::MUSTSET, "", ""))
    {
      readConfig(config, "agencyCode",       agencyCode,       Config::MUSTSET,  "TUG", "Identify the agency providing the data");
      readConfig(config, "timeStart",        timeStart,        Config::MUSTSET,  "", "Start time of the data");
      readConfig(config, "timeEnd",          timeEnd,          Config::MUSTSET,  "", "End time of the data ");
      readConfig(config, "observationCode",  observationCode,  Config::MUSTSET,  "P", "Technique used to generate the SINEX solution");
      readConfig(config, "solutionContents", solutionContents, Config::MUSTSET,  "", "Marker name if single station, MIX if multiple stations");
      readConfig(config, "description",      description,      Config::MUSTSET,  "", "Organizitions gathering/alerting the file contents");
      readConfig(config, "output",           output,           Config::MUSTSET,  "", "Description of the file contents");
      readConfig(config, "contact",          contact,          Config::MUSTSET,  "", "Address of the relevant contact e-mail");
      readConfig(config, "software",         software,         Config::MUSTSET,  "GROOPS", "Software used to generate the file");
      readConfig(config, "hardware",         hardware,         Config::MUSTSET,  "", "Computer hardware on which above software was run");
      readConfig(config, "input",            input,            Config::MUSTSET,  "", "Brief description of the input used to generate this solution");
      readConfig(config, "versionNumber",    versionNumber,    Config::MUSTSET,  "", "Unique identifier of the product, same as in file name, e.g. 000");
      readConfig(config, "inputfileComment", fileNameComment,  Config::OPTIONAL, "", "comments in the comment block from a file (truncated at 80 characters per line)");
      readConfig(config, "comment",          comment,          Config::OPTIONAL, "", "comments in the comment block");
      endSequence(config);
    }
    if(isCreateSchema(config)) return;

    // header line
    Time timeCurrent = System::now();

    // Reference time is mean value from start and end time
    timeRef = (timeStartObs + timeEndObs) * 0.5;

    // ==================================================

    // read data from files
    // --------------------

    // station informations
    logStatus << "reading station list from <" << fileNameStationList << ">" << Log::endl;
    std::vector<std::string> stationList;
    readFileStringList(fileNameStationList, stationList);

    logStatus << "reading station antenna definitions from <" << fileNameAntennaDef << ">" << Log::endl;
    std::vector<GnssAntennaDefinitionPtr> antennaDefinitionList;
    readFileGnssAntennaDefinition(fileNameAntennaDef, antennaDefinitionList);

    std::vector <Double> geoidHeightVector;
    if(!fileNameGeoidHeight.empty())
    {
      logStatus << "reading geoid heights <" << fileNameGeoidHeight << ">" << Log::endl;
      GriddedData griddedData;
      readFileGriddedData(fileNameGeoidHeight, griddedData);

      for(UInt i=0; i<stationList.size(); i++)
      {
        geoidHeightVector.push_back(griddedData.values.at(0).at(i));
      }
    }

    logStatus << "reading grid position <" << fileNameGridPos << ">" << Log::endl;
    GriddedData gridPos;
    readFileGriddedData(fileNameGridPos, gridPos);

    logStatus << "reading station infos from <" << fileNameStationInfo << ">" << Log::endl;
    std::vector<GnssStationInfo> stationInfos;
    VariableList fileNameVariableList;
    addVariable(variableStationName, fileNameVariableList);
    for(const auto &station : stationList)
    {
      fileNameVariableList[variableStationName]->setValue(station);
      GnssStationInfo stationInfo;
      readFileGnssStationInfo(fileNameStationInfo(fileNameVariableList), stationInfo);

      stationInfo.fillAntennaPattern(antennaDefinitionList);
      stationInfos.push_back(stationInfo);
    }

    // troposhere data
    std::vector<Matrix> tropDataStations(stationList.size()), tropSigmasStations(fileNameTropoSigmas.empty() ? 0 : stationList.size());
    if(!fileNameTropoData.empty())
    {
      logStatus << "reading troposphere data from <" << fileNameTropoData << ">" << Log::endl;

      for(UInt i = 0; i < stationList.size(); i++)
      {
        fileNameVariableList[variableStationName]->setValue(stationList.at(i));
        readFileMatrix(fileNameTropoData(fileNameVariableList),   tropDataStations.at(i));
        if(!fileNameTropoSigmas.empty())
          readFileMatrix(fileNameTropoSigmas(fileNameVariableList), tropSigmasStations.at(i));
      }
    }

    // ==================================================

    // write SINEX_TRO file
    // ----------------
    logStatus << "write SINEX_TRO file <" << fileNameTropoSinex << ">" << Log::endl;
    OutFile file(fileNameTropoSinex);

    // SINEX_TRO header line
    file << "%=TRO 2.00 " << std::setw(3) << agencyCode.substr(0, 3) << " " << time2str(timeCurrent) << " " << std::setw(3) << agencyCode.substr(0, 3);
    file << " " << time2str(timeStart) << " " << time2str(timeEnd) << " " << std::setw(1) << observationCode.substr(0, 1);
    file << " " << solutionContents.substr(0, 4) << std::endl;

    // Block: FILE/REFERENCE
    file << "*-------------------------------------------------------------------------------" << std::endl;
    file << "+FILE/REFERENCE" << std::endl;
    file << "*INFO_TYPE_________ INFO________________________________________________________" << std::endl;
    file << " DESCRIPTION        " << resize(description, 60) << std::endl;
    file << " OUTPUT             " << resize(output, 60) << std::endl;
    file << " CONTACT            " << resize(contact, 60) << std::endl;
    file << " SOFTWARE           " << resize(software, 60) << std::endl;
    file << " HARDWARE           " << resize(hardware, 60) << std::endl;
    file << " INPUT              " << resize(input, 60) << std::endl;
    file << " VERSION NUMBER     " << resize(versionNumber, 60) << std::endl;
    file << "-FILE/REFERENCE" << std::endl;

    // Block: TROP/DESCRIPTION
    std::vector<std::string> columns = {"TROTOT", "TGNTOT", "TGETOT"};
    if(!fileNameTropoSigmas.empty())
      columns = {"TROTOT", "STDDEV", "TGNTOT", "STDDEV", "TGETOT", "STDDEV"};
    file << "*-------------------------------------------------------------------------------" << std::endl;
    file << "+TROP/DESCRIPTION" << std::endl;
    file << "*_________KEYWORD____________ ___VALUES(S)______________________________________" << std::endl;
    file << " TROPO PARAMETER NAMES       "; for(const auto &c : columns) {file << " " << c;}; file << std::endl;
    file << " TROPO PARAMETER UNITS       "; for(UInt i = 0; i < columns.size(); i++) {file << "  1e+03";}; file << std::endl;
    file << " TROPO PARAMETER WIDTH       "; for(UInt i = 0; i < columns.size(); i++) {file << "      6";}; file << std::endl;
    file << " TROPO MODELING METHOD        " << tropoModelingMethod << std::endl;
    file << " TROPO SAMPLING INTERVAL      " << troSampling%"%f"s << std::endl;
    if(!metDataSource.empty())
      file << " SOURCE OF MET/DATA           " << metDataSource << std::endl;
    file << " A PRIORI TROPOSPHERE         " << aPrioriTropoModel << std::endl;
    file << " TROPO MAPPING FUNCTION       " << tropoMappingFunction << std::endl;
    file << " GRADS MAPPING FUNCTION       " << gradientMappingFunction << std::endl;
    file << " DATA SAMPLING INTERVAL       " << dataSampling%"%f"s << std::endl;
    file << " ELEVATION CUTOFF ANGLE       " << elevationCutoff%"%f"s << std::endl;
    file << " OBSERVATION WEIGHTING        " << observationWeighting << std::endl;
    file << " GNSS SYSTEMS                 " << gnssSystems << std::endl;
    file << " TIME SYSTEM                  " << timeSystem << std::endl;
    if(!oceanTideModel.empty())
      file << " OCEAN TIDE LOADING MODEL     " << oceanTideModel << std::endl;
    if(!atmosphericTideModel.empty())
      file << " ATMOSPH TIDE LOADING MODEL   " << atmosphericTideModel << std::endl;
    if(!geoidModel.empty())
      file << " GEOID MODEL                  " << geoidModel << std::endl;
    file << "-TROP/DESCRIPTION" << std::endl;

    // Block: FILE/COMMENT
    // comment block
    if(comment.size() || !fileNameComment.empty())
    {
      file << "*-------------------------------------------------------------------------------" << std::endl;
      file << "+FILE/COMMENT" << std::endl;
      if(!fileNameComment.empty())
      {
        InFile commentFile(fileNameComment);
        std::string line;
        while(std::getline(commentFile, line))
          file << " " << line << std::endl;
      }
      for(const auto &line : comment)
        file << " " << line << std::endl;
      file << "-FILE/COMMENT" << std::endl;
    }

    // Block: SITE/ID
    blockText = SiteId(stationInfos, gridPos, geoidHeightVector);
    file << "*-------------------------------------------------------------------------------" << std::endl;
    file << "+SITE/ID" << std::endl;
    file << "*STATION__ PT __DOMES__ T _STATION_DESCRIPTION__ _LONGITUDE _LATITUDE_ _HGT_ELI_ _HGT_MSL_" << std::endl;
    file << blockText;
    file << "-SITE/ID" << std::endl;

    // Block: SITE/COORDINATES
    blockText = SiteCoordinates(stationInfos, gridPos, timeStart, timeEnd, systemCode, remark);
    file << "*-------------------------------------------------------------------------------" << std::endl;
    file << "+SITE/COORDINATES" << std::endl;
    file << "*STATION__ PT SOLN T __DATA_START__ __DATA_END____ __STA_X_____  __STA_Y_____  __STA_Z_____ SYSTEM REMRK" << std::endl;
    file << blockText;
    file << "-SITE/COORDINATES" << std::endl;

    // BLOCK: SITE/RECEIVER
    blockText = SiteReceiver(stationInfos, timeRef);
    file << "*-------------------------------------------------------------------------------" << std::endl;
    file << "+SITE/RECEIVER" << std::endl;
    file << "*STATION__ PT SOLN T __DATA_START__ __DATA_END____ DESCRIPTION_________ S/N_________________ FIRMW______" << std::endl;
    file << blockText;
    file << "-SITE/RECEIVER" << std::endl;

    // BLOCK SITE/ANTENNA
    blockText = SiteAntenna(stationInfos, timeRef, antennaModel);
    file << "*-------------------------------------------------------------------------------" << std::endl;
    file << "+SITE/ANTENNA" << std::endl;
    file << "*STATION__ PT SOLN T __DATA_START__ __DATA_END____ DESCRIPTION_________ S/N_________________ PCV_MODEL_" << std::endl;
    file << blockText;
    file << "-SITE/ANTENNA" << std::endl;

    // BLOCK SITE/ECCENTRICITY
    blockText = SiteEccentricity(stationInfos, timeRef);
    file << "*-------------------------------------------------------------------------------" << std::endl;
    file << "+SITE/ECCENTRICITY" << std::endl;
    file << "*" << std::setw(80) << "UP______ NORTH___ EAST____" << std::endl;
    file << "*STATION__ PT SOLN T __DATA_START__ __DATA_END____ AXE ARP->BENCHMARK(M)_________" << std::endl;
    file << blockText;
    file << "-SITE/ECCENTRICITY" << std::endl;

    // BLOCK: TROP/SOLUTION
    blockText = TropSolution(stationInfos, tropDataStations, tropSigmasStations);
    file << "*-------------------------------------------------------------------------------" << std::endl;
    file << "+TROP/SOLUTION" << std::endl;
    file << "*STATION__ ____EPOCH_____"; for(const auto &c : columns) {file << " " << c;}; file << std::endl;
    file << blockText;
    file << "-TROP/SOLUTION" << std::endl;

    file << "%=ENDTRO" << std::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string GnssTroposphere2TropoSinex::SiteId(const std::vector<GnssStationInfo> &stationInfos, const GriddedData &griddedData, const std::vector<Double> &geoidHeightVector)
{
  try
  {
    std::stringstream ss;
    Ellipsoid ellipsoid;

    UInt i = 0;
    for(const auto &stationInfo : stationInfos)
    {
      Angle lon, lat;
      Double ellipsoidHeight;
      ellipsoid(griddedData.points.at(i), lon, lat, ellipsoidHeight);

      ss << " " << String::upperCase(resize(stationInfo.markerName, 9)) << "  A " << resize(stationInfo.markerNumber, 9) << " P " << resize(stationInfo.comment, 22) << " "
         << (std::fmod(lon+2*PI, 2*PI)*RAD2DEG)%"%10.6f "s << (lat*RAD2DEG)%"%10.6f "s << ellipsoidHeight%"%9.3f "s << (geoidHeightVector.empty() ? "" : geoidHeightVector.at(i) % "%9.3f"s) << std::endl;

      i++;
     }

     return ss.str();
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string GnssTroposphere2TropoSinex::SiteCoordinates(const std::vector<GnssStationInfo> &stationInfos, const GriddedData &griddedData, const Time &timeStart,
                                                        const Time &timeEnd, std::string &systemCode, std::string &remark)
{
  try
  {
    std::stringstream ss;

    UInt i = 0;
    for(const auto &stationInfo : stationInfos)
    {
      ss << " " << String::upperCase(resize(stationInfo.markerName, 9)) << "  A    1 P " << time2str(timeStart) << " " << time2str(timeEnd) << " " << griddedData.points.at(i).x() % "%12.3f "s;
      ss << " " << griddedData.points.at(i).y()%"%12.3f "s << " " << griddedData.points.at(i).z()%"%12.3f "s << resize(systemCode, 6) << " " << resize(remark, 5) << std::endl;

      i++;
    }

    return ss.str();
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string GnssTroposphere2TropoSinex::SiteReceiver(const std::vector<GnssStationInfo> &stationInfos, const Time &timeRef)
{
  try
  {
    std::stringstream ss;

    for(const auto &stationInfo : stationInfos)
    {
      const UInt idRecv = stationInfo.findReceiver(timeRef);
      if(idRecv == NULLINDEX)
      {
        logWarning << stationInfo.markerName << ": no receiver found at " << timeRef.dateTimeStr() << Log::endl;
        continue;
      }
      const GnssReceiverInfo recv = stationInfo.receiver.at(idRecv);
      ss << " " << String::upperCase(resize(stationInfo.markerName, 9)) << "  A    1 P " << time2str(recv.timeStart) << " " << time2str(recv.timeEnd) << " " << resize(recv.name, 20)
         << " " << resize(recv.serial.empty() ? std::string(20, '-') : recv.serial, 20) << " " << resize(recv.version.empty() ? std::string(11, '-') : recv.version, 11) << std::endl;
    }

    return ss.str();
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string GnssTroposphere2TropoSinex::SiteAntenna(const std::vector<GnssStationInfo> &stationInfos, const Time &timeRef, std::string &antennaModel)
{
  try
  {
    std::stringstream ss;

    for(const auto &stationInfo : stationInfos)
    {
      const UInt idAnt = stationInfo.findAntenna(timeRef);
      if(idAnt == NULLINDEX)
      {
        logWarning << stationInfo.markerName << ": no antenna found at " << timeRef.dateTimeStr() << Log::endl;
        continue;
      }
      const GnssAntennaInfo ant = stationInfo.antenna.at(idAnt);

      ss << " " << String::upperCase(resize(stationInfo.markerName, 9)) << "  A    1 P " << time2str(ant.timeStart) << " " << time2str(ant.timeEnd) << " "
         << resize(ant.name, 15) << " " << (ant.radome.empty() ? "NONE" : resize(ant.radome, 4)) << " " << resize(ant.serial.empty() ? std::string(20, '-') : ant.serial, 20) << " "
         << resize(antennaModel, 10) << std::endl;
    }

    return ss.str();
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string GnssTroposphere2TropoSinex::SiteEccentricity(const std::vector<GnssStationInfo> &stationInfos, const Time &timeRef)
{
  try
  {
    std::stringstream ss;

    for(const auto &stationInfo : stationInfos)
    {
      const UInt idAnt = stationInfo.findAntenna(timeRef);
      if(idAnt == NULLINDEX)
      {
        logWarning << stationInfo.markerName << ": no antenna found at " << timeRef.dateTimeStr() << Log::endl;
        continue;
      }
      const GnssAntennaInfo ant = stationInfo.antenna.at(idAnt);

      ss <<" " << String::upperCase(resize(stationInfo.markerName, 9)) << "  A    1 P " << time2str(ant.timeStart) << " "
         << time2str(ant.timeEnd) << " UNE " << ant.position.z() % "%8.4f "s << ant.position.x() % "%8.4f "s << ant.position.y() % "%8.4f"s << std::endl;
    }

    return ss.str();
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string GnssTroposphere2TropoSinex::TropSolution(const std::vector<GnssStationInfo> &stationInfos, const std::vector<Matrix> &tropDataMatrix, const std::vector<Matrix> &tropSigmasMatrix)
{
  try
  {
    std::stringstream ss;

    for(UInt idStation = 0; idStation < stationInfos.size(); idStation++)
      for(UInt idEpoch = 0; idEpoch < tropDataMatrix.at(idStation).rows(); idEpoch++)
      {
        Vector epochData = 1e3 * tropDataMatrix.at(idStation).row(idEpoch).trans();
        ss << " " << String::upperCase(resize(stationInfos.at(idStation).markerName, 9)) << " " << time2str(mjd2time(epochData(0)*1e-3));
        if(tropSigmasMatrix.empty())
          ss << (epochData(1)+epochData(2))%" %6.1f"s << (epochData(3)+epochData(4))%" %6.2f"s << (epochData(5)+epochData(6))%" %6.2f"s << std::endl;
        else
        {
          Vector epochSigmas = 1e3 * tropSigmasMatrix.at(idStation).row(idEpoch).trans();
          ss << (epochData(1)+epochData(2))%" %6.1f"s << epochSigmas(1)%" %6.1f"s
             << (epochData(3)+epochData(4))%" %6.2f"s << epochSigmas(2)%" %6.2f"s
             << (epochData(5)+epochData(6))%" %6.2f"s << epochSigmas(3)%" %6.2f"s << std::endl;
        }
      }

    return ss.str();
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
