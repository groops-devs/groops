/***********************************************/
/**
* @file gnssNormals2Sinex.cpp
*
* @brief Write GNSS data/metadata and normal equations to SINEX format.
*
* @author Sebastian Strasser
* @date 2019-05-21
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Write GNSS data/metadata and \file{normal equations}{normalEquation} to
\href{http://www.iers.org/IERS/EN/Organization/AnalysisCoordinator/SinexFormat/sinex.html}{SINEX format}.

Normal equations usually come from \program{GnssProcessing}
(e.g. from \reference{GNSS satellite orbit determination and station network analysis}{cookbook.gnssNetwork}).
Metadata input files include \configFile{stationInfo/transmitterInfo}{platform}, \configFile{antennaDefinition}{gnssAntennaDefinition},
and \configFile{stationList/transmitterList}{stringList}, see \program{GnssAntex2AntennaDefinition}.

See also \program{Sinex2Normals} and \program{NormalsSphericalHarmonics2Sinex}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "inputOutput/fileSinex.h"
#include "files/filePlatform.h"
#include "files/fileMatrix.h"
#include "files/fileNormalEquation.h"
#include "files/fileStringTable.h"

/***** CLASS ***********************************/

/** @brief Write GNSS data/metadata and normal equations to SINEX format.
* @ingroup programsConversionGroup */
class GnssNormals2Sinex
{
public:
  class TransmitterConstellation
  {
  public:
    FileName fileNameTransmitterList, fileNameTransmitterInfo, fileNameAntennaDef;
    std::string variablePrn;
  };

  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GnssNormals2Sinex, SINGLEPROCESS, "Write GNSS data/metadata and normal equations to SINEX format.", Conversion, Gnss, NormalEquation)

/***********************************************/

template<> Bool readConfig(Config &config, const std::string &name, GnssNormals2Sinex::TransmitterConstellation &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  if(!readConfigSequence(config, name, mustSet, defaultValue, annotation))
    return FALSE;
  readConfig(config, "inputfileTransmitterList",   var.fileNameTransmitterList, Config::MUSTSET,  "",    "transmitter PRNs used in solution");
  readConfig(config, "inputfileTransmitterInfo",   var.fileNameTransmitterInfo, Config::MUSTSET,  "{groopsDataDir}/gnss/transmitter/transmitterInfo/igs/igs20/transmitterInfo_igs20.{prn}.xml", "transmitter info file template");
  readConfig(config, "inputfileAntennaDefinition", var.fileNameAntennaDef,      Config::MUSTSET,  "{groopsDataDir}/gnss/transmitter/antennaDefinition/igs/igs20/transmitterInfo_igs20.dat",     "transmitter phase centers and variations (ANTEX)");
  readConfig(config, "variablePrn",                var.variablePrn,             Config::DEFAULT,  "prn", "loop variable for PRNs from transmitter list");
  endSequence(config);
  return TRUE;
}

/***********************************************/

void GnssNormals2Sinex::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName    fileNameSinexNormals, fileNameSinexCoords, fileNameNormals, fileNameSolution, fileNameSigmax, fileNameApriori, fileNameAprioriSigma, fileNameAprMat;
    FileName    fileNameStationList, fileNameStationInfo, fileNameAntennaDef;
    Time        timeRef, timeStartObs, timeEndObs;
    std::string variableStationName, antennaModel;
    Double      sampling = 0;
    std::vector<TransmitterConstellation> constellations;
    Sinex       sinex;

    readConfig(config, "outputfileSinexNormals",     fileNameSinexNormals, Config::OPTIONAL, "", "full SINEX file including normal equations");
    readConfig(config, "outputfileSinexCoordinates", fileNameSinexCoords,  Config::OPTIONAL, "", "SINEX file without normal equations (station coordinates file)");
    readConfig(config, "inputfileNormals",           fileNameNormals,      Config::MUSTSET,  "", "normal equation matrix");
    readConfig(config, "inputfileSolution",          fileNameSolution,     Config::OPTIONAL, "", "parameter vector");
    readConfig(config, "inputfileSigmax",            fileNameSigmax,       Config::OPTIONAL, "", "standard deviations of the parameters (sqrt of the diagonal of the inverse normal equation)");
    readConfig(config, "inputfileApriori",           fileNameApriori,      Config::MUSTSET,  "", "apriori parameter vector");
    readConfig(config, "inputfileAprioriSigma",      fileNameAprioriSigma, Config::OPTIONAL, "", "constraint sigmas for apriori parameter vector");
    readConfig(config, "inputfileAprioriMatrix",     fileNameAprMat,       Config::OPTIONAL, "", "normal equation matrix of applied constraints");
    readConfig(config, "transmitterConstellation",   constellations,       Config::MUSTSET,  "", "transmitter constellation metadata");
    if(readConfigSequence(config, "stations", Config::MUSTSET, "", ""))
    {
      readConfig(config, "inputfileStationList",       fileNameStationList, Config::MUSTSET,  "", "stations contained in normal equations");
      readConfig(config, "inputfileStationInfo",       fileNameStationInfo, Config::MUSTSET,  "{groopsDataDir}/gnss/receiverStation/stationInfo/igs/stationInfo.{station}.xml", "station info file template");
      readConfig(config, "inputfileAntennaDefinition", fileNameAntennaDef,  Config::MUSTSET,  "{groopsDataDir}/gnss/receiverStation/antennaDefinition/igs/igs20/antennaDefinition_igs20.dat", "station phase centers and variations (ANTEX)");
      readConfig(config, "variableStationName",        variableStationName, Config::DEFAULT,  "station", "loop variable for station names from station list");
      readConfig(config, "observationTimeStart",       timeStartObs,        Config::MUSTSET,  "", "start time for which solution has observations");
      readConfig(config, "observationTimeEnd",         timeEndObs,          Config::MUSTSET,  "", "end time for which solution has observations");
      endSequence(config);
    }
    readConfig(config, "time",                       timeRef,             Config::MUSTSET,  "", "reference time for parameters");
    readConfig(config, "sampling",                   sampling,            Config::OPTIONAL, "", "[seconds] observation sampling");
    readConfig(config, "antennaCalibrationModel",    antennaModel,        Config::MUSTSET,  "", "e.g. IGS14_WWWW (WWWW = ANTEX release GPS week)");
    readConfig(config, "sinexHeader",                sinex,               Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    // ==================================================

    // read data from files
    // --------------------
    Vector x;
    if(!fileNameSolution.empty())
    {
      logStatus<<"reading solution from <"<<fileNameSolution<<">"<<Log::endl;
      readFileMatrix(fileNameSolution, x);
    }

    Vector sigmax;
    if(!fileNameSigmax.empty())
    {
      logStatus<<"reading standard deviations from <"<<fileNameSigmax<<">"<<Log::endl;
      readFileMatrix(fileNameSigmax, sigmax);
    }

    Vector x0;
    if(!fileNameApriori.empty())
    {
      logStatus<<"reading apriori solution from <"<<fileNameApriori<<">"<<Log::endl;
      readFileMatrix(fileNameApriori, x0);
    }

    Vector sigmax0;
    if(!fileNameAprioriSigma.empty())
    {
      logStatus<<"reading constraint sigmas from <"<<fileNameAprioriSigma<<">"<<Log::endl;
      readFileMatrix(fileNameAprioriSigma, sigmax0);
    }

    Matrix N, n;
    NormalEquationInfo info;
    logStatus<<"reading normal equation matrix from <"<<fileNameNormals<<">"<<Log::endl;
    readFileNormalEquation(fileNameNormals, info, N, n);
    fillSymmetric(N);
    const UInt countParameter = N.rows();

    Matrix dN;
    std::vector<Bool> parameterIsConstrained(countParameter, FALSE);
    if(!fileNameAprMat.empty())
    {
      logStatus<<"reading normal equation matrix of applied constraints <"<<fileNameAprMat<<">"<<Log::endl;
      Vector n;
      NormalEquationInfo info;
      readFileNormalEquation(fileNameAprMat, info, dN, n);
      fillSymmetric(dN);
      if(dN.rows() != parameterIsConstrained.size())
        throw(Exception("Parameter count in constraint matrix and normal equation matrix differs (" + dN.rows()%"%i"s + " vs. "+ N.rows()%"%i"s +" )."));
      for(UInt i=0; i<dN.rows(); i++)
        parameterIsConstrained.at(i) = (dN(i, i) != 0.0);
    }

    logStatus<<"reading station list from <"<<fileNameStationList<<">"<<Log::endl;
    std::vector<std::string> stationList;
    readFileStringList(fileNameStationList, stationList);

    logStatus<<"reading station antenna definitions from <"<<fileNameAntennaDef<<">"<<Log::endl;
    std::vector<GnssAntennaDefinitionPtr> antennaDefinitionList;
    readFileGnssAntennaDefinition(fileNameAntennaDef, antennaDefinitionList);

    logStatus<<"reading station infos from <"<<fileNameStationInfo<<">"<<Log::endl;
    std::vector<Platform> stationInfos;
    VariableList fileNameVariableList;
    for(const auto &station : stationList)
    {
      fileNameVariableList.setVariable(variableStationName, station);
      Platform stationInfo;
      readFilePlatform(fileNameStationInfo(fileNameVariableList), stationInfo);
      stationInfo.fillGnssAntennaDefinition(antennaDefinitionList);
      stationInfos.push_back(stationInfo);
    }

    std::vector<Platform> transmitterInfos;
    for(const auto &constellation : constellations)
    {
      logStatus<<"reading transmitter list from <"<<constellation.fileNameTransmitterList<<">"<<Log::endl;
      std::vector<std::string> transmitterList;
      readFileStringList(constellation.fileNameTransmitterList, transmitterList);

      logStatus<<"reading transmitter antenna definitions from <"<<constellation.fileNameAntennaDef<<">"<<Log::endl;
      std::vector<GnssAntennaDefinitionPtr> antennaDefinitionList;
      readFileGnssAntennaDefinition(constellation.fileNameAntennaDef, antennaDefinitionList);

      logStatus<<"reading transmitter infos from <"<<constellation.fileNameTransmitterInfo<<">"<<Log::endl;
      for(const auto &prn : transmitterList)
      {
        fileNameVariableList.setVariable(constellation.variablePrn, prn);
        Platform transmitterInfo;
        readFilePlatform(constellation.fileNameTransmitterInfo(fileNameVariableList), transmitterInfo);
        transmitterInfo.fillGnssAntennaDefinition(antennaDefinitionList);
        transmitterInfos.push_back(transmitterInfo);
      }
    }

    // ==================================================

    // add data to SINEX
    // -----------------
    {
      SinexBlockPtr block = sinex.addBlock("SITE/ID");
      *block<<"*SITE PT __DOMES__ T _STATION DESCRIPTION__ _LONGITUDE_ _LATITUDE__ HEIGHT_"<<std::endl;

      auto rad2DegMinSec = [] (Double rad, Double &deg, Double &min, Double& sec)
      {
        deg = std::floor(rad*RAD2DEG);
        min = std::floor((rad*RAD2DEG-deg)*60);
        sec = ((rad*RAD2DEG-deg)*60-min)*60;
      };

      Ellipsoid ellipsoid;
      for(const auto &stationInfo : stationInfos)
      {
        Angle lon, lat;
        Double height, lonDeg, lonMin, lonSec, latDeg, latMin, latSec;
        ellipsoid(stationInfo.approxPosition, lon, lat, height);
        rad2DegMinSec(lon < 0 ? lon+2*PI : lon, lonDeg, lonMin, lonSec);
        rad2DegMinSec(lat, latDeg, latMin, latSec);

        std::stringstream ss;
        *block<<String::upperCase(Sinex::resize(stationInfo.markerName, 4))<<"  A "<<Sinex::resize(stationInfo.markerNumber, 9)<<" P "<<Sinex::resize(stationInfo.comment, 22)<<" "
              <<lonDeg%"%3i "s<<lonMin%"%2i "s<<lonSec%"%4.1f "s<<latDeg%"%3i "s<<latMin%"%2i "s<<latSec%"%4.1f "s<<height%"%7.1f"s<<std::endl;
      }
    }
    // -----------------
    {
      SinexBlockPtr block = sinex.addBlock("SITE/ECCENTRICITY");
      *block<<"*SITE PT SOLN T _DATA START_ __DATA_END__ AXE _ECC_U__ _ECC_N__ _ECC_E__"<<std::endl;
      for(const auto &stationInfo : stationInfos)
      {
        auto ant = stationInfo.findEquipment<PlatformGnssAntenna>(timeRef);
        if(!ant)
        {
          logWarning<<stationInfo.markerName<<": no antenna found at "<<timeRef.dateTimeStr()<<Log::endl;
          continue;
        }
        *block<<String::upperCase(Sinex::resize(stationInfo.markerName, 4))<<"  A    1 P "<<Sinex::time2str(ant->timeStart)<<" "
              <<Sinex::time2str(ant->timeEnd)<<" UNE "<<ant->position.z()%"%8.4f "s<<ant->position.x()%"%8.4f "s<<ant->position.y()%"%8.4f "s<<std::endl;
      }
    }
    // -----------------
    {
      SinexBlockPtr block = sinex.addBlock("SITE/RECEIVER");
      *block<<"*SITE PT SOLN T _DATA START_ __DATA_END__ ___RECEIVER_TYPE____ _S/N_ _FIRMWARE__"<<std::endl;
      for(const auto &stationInfo : stationInfos)
      {
        auto recv = stationInfo.findEquipment<PlatformGnssReceiver>(timeRef);
        if(!recv)
        {
          logWarning<<stationInfo.markerName<<": no receiver found at "<<timeRef.dateTimeStr()<<Log::endl;
          continue;
        }
        *block<<String::upperCase(Sinex::resize(stationInfo.markerName, 4))<<"  A    1 P "<<Sinex::time2str(recv->timeStart)<<" "<<Sinex::time2str(recv->timeEnd)<<" "
              <<Sinex::resize(recv->name, 20)<<" "<<(recv->serial.empty() ? "-----" : Sinex::resize(recv->serial, 5))
              <<" "<<(recv->version.empty() ? std::string(11, '-') : Sinex::resize(recv->version, 11))<<std::endl;
      }
    }
    // -----------------
    {
      SinexBlockPtr block = sinex.addBlock("SITE/ANTENNA");
      *block<<"*SITE PT SOLN T _DATA START_ __DATA_END__ ____ANTENNA_TYPE____ _S/N_"<<std::endl;
      for(const auto &stationInfo : stationInfos)
      {
        auto ant = stationInfo.findEquipment<PlatformGnssAntenna>(timeRef);
        if(!ant)
        {
          logWarning<<stationInfo.markerName<<": no antenna found at "<<timeRef.dateTimeStr()<<Log::endl;
          continue;
        }
        Angle roll, pitch, yaw;
        Rotary3d(ant->local2antennaFrame.matrix()).cardan(roll, pitch, yaw);
        *block<<String::upperCase(Sinex::resize(stationInfo.markerName, 4))<<"  A    1 P "<<Sinex::time2str(ant->timeStart)<<" "<<Sinex::time2str(ant->timeEnd)<<" "
              <<Sinex::resize(ant->name, 15)<<" "<<(ant->radome.empty() ? "NONE" : Sinex::resize(ant->radome, 4))<<" "<<(ant->serial.empty() ? "-----" : Sinex::resize(ant->serial, 5))
              <<" "<<(yaw*RAD2DEG)%"% 4i"s<<std::endl;
      }
    }
    // -----------------
    {
      SinexBlockPtr block = sinex.addBlock("SITE/GPS_PHASE_CENTER");
      *block<<"*____ANTENNA_TYPE____ _S/N_ _L1_U_ _L1_N_ _L1_E_ _L2_U_ _L2_N_ _L2_E_ _ANTMODEL_"<<std::endl;
      for(const auto &stationInfo : stationInfos)
      {
        auto ant = stationInfo.findEquipment<PlatformGnssAntenna>(timeRef);
        if(!ant)
        {
          logWarning<<stationInfo.markerName<<": no antenna found at "<<timeRef.dateTimeStr()<<Log::endl;
          continue;
        }
        auto p1 = std::find_if(ant->antennaDef->patterns.begin(), ant->antennaDef->patterns.end(), [](const GnssAntennaPattern &pattern) {return pattern.type == GnssType::L1_G;});
        auto p2 = std::find_if(ant->antennaDef->patterns.begin(), ant->antennaDef->patterns.end(), [](const GnssAntennaPattern &pattern) {return pattern.type == GnssType::L2_G;});
        if(p1 == ant->antennaDef->patterns.end() || p2 == ant->antennaDef->patterns.end())
        {
          logWarning<<stationInfo.markerName <<": GPS phase center not found for antenna "<<ant->name<<" "<<ant->radome<<" "<<ant->serial<<Log::endl;
          continue;
        }
        *block<<Sinex::resize(ant->name, 15)<<" "<<(ant->radome.empty() ? "NONE" : Sinex::resize(ant->radome, 4))<<" "<<(ant->serial.empty() ? "-----" : Sinex::resize(ant->serial, 5))<<" "
              <<Sinex::format(p1->offset.z())<<" "<<Sinex::format(p1->offset.x())<<" "<<Sinex::format(p1->offset.y())<<" "
              <<Sinex::format(p2->offset.z())<<" "<<Sinex::format(p2->offset.x())<<" "<<Sinex::format(p2->offset.y())<<" "
              <<Sinex::resize(antennaModel, 10)<<std::endl;
      }
    }
    // -----------------
    {
      SinexBlockPtr block = sinex.addBlock("SITE/GAL_PHASE_CENTER");
      *block<<"*____ANTENNA_TYPE____ _S/N_ _L1_U_ _L1_N_ _L1_E_ _L5_U_ _L5_N_ _L5_E_ _ANTMODEL_"<<std::endl;
      *block<<"*____ANTENNA_TYPE____ _S/N_ _L6_U_ _L6_N_ _L6_E_ _L7_U_ _L7_N_ _L7_E_ _ANTMODEL_"<<std::endl;
      *block<<"*____ANTENNA_TYPE____ _S/N_ _L8_U_ _L8_N_ _L8_E_ ____________________ _ANTMODEL_"<<std::endl;
      for(const auto &stationInfo : stationInfos)
      {
        auto ant = stationInfo.findEquipment<PlatformGnssAntenna>(timeRef);
        if(!ant)
        {
          logWarning<<stationInfo.markerName<<": no antenna found at "<<timeRef.dateTimeStr()<<Log::endl;
          continue;
        }
        const auto &patterns = ant->antennaDef->patterns;
        auto p1 = std::find_if(patterns.begin(), patterns.end(), [](const GnssAntennaPattern &pattern){return pattern.type == GnssType::L1_E;});
        auto p5 = std::find_if(patterns.begin(), patterns.end(), [](const GnssAntennaPattern &pattern){return pattern.type == GnssType::L5_E;});
        auto p6 = std::find_if(patterns.begin(), patterns.end(), [](const GnssAntennaPattern &pattern){return pattern.type == GnssType::L6_E;});
        auto p7 = std::find_if(patterns.begin(), patterns.end(), [](const GnssAntennaPattern &pattern){return pattern.type == GnssType::L7_E;});
        auto p8 = std::find_if(patterns.begin(), patterns.end(), [](const GnssAntennaPattern &pattern){return pattern.type == GnssType::L8_E;});
        if(p1 == patterns.end() || p5 == patterns.end() || p6 == patterns.end() || p7 == patterns.end() || p8 == patterns.end())
          continue;
        *block<<Sinex::resize(ant->name, 15)<<" "<<(ant->radome.empty() ? "NONE" : Sinex::resize(ant->radome, 4))<<" "<<(ant->serial.empty() ? "-----" : Sinex::resize(ant->serial, 5))<<" "
              <<Sinex::format(p1->offset.z())<<" "<<Sinex::format(p1->offset.x())<<" "<<Sinex::format(p1->offset.y())<<" "
              <<Sinex::format(p5->offset.z())<<" "<<Sinex::format(p5->offset.x())<<" "<<Sinex::format(p5->offset.y())<<" "<<Sinex::resize(antennaModel, 10)<<std::endl;
        *block<<Sinex::resize(ant->name, 15)<<" "<<(ant->radome.empty() ? "NONE" : Sinex::resize(ant->radome, 4))<<" "<<(ant->serial.empty() ? "-----" : Sinex::resize(ant->serial, 5))<<" "
              <<Sinex::format(p6->offset.z())<<" "<<Sinex::format(p6->offset.x())<<" "<<Sinex::format(p6->offset.y())<<" "
              <<Sinex::format(p7->offset.z())<<" "<<Sinex::format(p7->offset.x())<<" "<<Sinex::format(p7->offset.y())<<" "<<Sinex::resize(antennaModel, 10)<<std::endl;
        *block<<Sinex::resize(ant->name, 15)<<" "<<(ant->radome.empty() ? "NONE" : Sinex::resize(ant->radome, 4))<<" "<<(ant->serial.empty() ? "-----" : Sinex::resize(ant->serial, 5))<<" "
              <<Sinex::format(p8->offset.z())<<" "<<Sinex::format(p8->offset.x())<<" "<<Sinex::format(p8->offset.y())<<" "
              <<std::string(20, ' ')<<" "<<Sinex::resize(antennaModel, 10)<<std::endl;
      }
    }
    // -----------------
    {
      SinexBlockPtr block = sinex.addBlock("SATELLITE/ID");
      *block<<"*SVN_ PR COSPAR_ID T _DATA_START_ __DATA_END__ ______ANTENNA_______"<<std::endl;

      std::map<std::string, Time> svn2TimeStart, svn2TimeEnd;
      for(const auto &transmitterInfo : transmitterInfos)
        for(const auto &instrument : transmitterInfo.equipments)
        {
          auto ant = std::dynamic_pointer_cast<PlatformGnssAntenna>(instrument);
          if(ant)
          {
            svn2TimeStart[ant->serial] = (svn2TimeStart[ant->serial] == Time() ? ant->timeStart : std::min(svn2TimeStart[ant->serial], ant->timeStart));
            svn2TimeEnd[ant->serial]   = std::max(svn2TimeStart[ant->serial], ant->timeEnd);
          }
        }

      for(const auto &transmitterInfo : transmitterInfos)
      {
        auto ant = transmitterInfo.findEquipment<PlatformGnssAntenna>(timeRef);
        if(!ant)
        {
          logWarning<<transmitterInfo.markerName<<": no antenna found at "<<timeRef.dateTimeStr()<<Log::endl;
          continue;
        }
        *block<<ant->serial<<" "<<transmitterInfo.markerNumber.substr(1,2)<<" "<<ant->radome<<" P "<<Sinex::time2str(svn2TimeStart[ant->serial])
              <<" "<<Sinex::time2str(svn2TimeEnd[ant->serial])<<" "<<ant->name<<std::endl;
      }
    }
    // -----------------
    {
      SinexBlockPtr block = sinex.addBlock("SATELLITE/PHASE_CENTER");
      *block<<"*SVN_ L SATA_Z SATA_X SATA_Y L SATA_Z SATA_X SATA_Y _ANTMODEL_ T M"<<std::endl;

      auto antennaOffsetStr = [&](const PlatformGnssAntenna &ant, const GnssAntennaPattern &pattern)
      {
        Vector3d offset = ant.position + ant.local2antennaFrame.inverseTransform(pattern.offset);
        return pattern.type.str().substr(1,1)+" "+Sinex::format(offset.z())+" "+Sinex::format(offset.x())+" "+Sinex::format(offset.y());
      };

      for(const auto &transmitterInfo : transmitterInfos)
      {
        auto ant = transmitterInfo.findEquipment<PlatformGnssAntenna>(timeRef);
        if(!ant)
        {
          logWarning<<transmitterInfo.markerName<<": no antenna found at "<<timeRef.dateTimeStr()<<Log::endl;
          continue;
        }
        std::vector<const GnssAntennaPattern*> patterns;
        for(const auto &pattern : ant->antennaDef->patterns)
          if(pattern.type == GnssType::PHASE)
            patterns.push_back(&pattern);
        std::sort(patterns.begin(), patterns.end(), [](const auto &p1, const auto &p2){return p1->type < p2->type;});
        for(UInt i=0; i<patterns.size(); i+=2)
        {
          *block<<ant->serial<<" "<<antennaOffsetStr(*ant, *patterns.at(i))<<" "<<(i+1 < patterns.size() ? antennaOffsetStr(*ant, *patterns.at(i+1)) : std::string(22, ' '))
                <<" "<<Sinex::resize(antennaModel, 10)<<" A "<<(patterns.at(i)->pattern.rows() > 1 ? "F" : "E")<<std::endl;
        }
      }
    }
    // -----------------
    {
      SinexBlockPtr block = sinex.addBlock("SATELLITE/PHASE_CENTER");
      *block<<"*SITE PT SOLN T _DATA_START_ __DATA_END__ _MEAN_EPOCH_"<<std::endl;
      for(const auto &stationInfo : stationInfos)
        *block<<String::upperCase(Sinex::resize(stationInfo.markerName, 4))<<"  A    1 P "
              <<Sinex::time2str(timeStartObs)<<" "<<Sinex::time2str(timeEndObs)<<" "<<Sinex::time2str(0.5*(timeStartObs+timeEndObs))<<std::endl;
    }
    // -----------------
    {
      SinexBlockPtr block = sinex.addBlock("SOLUTION/STATISTICS");
      *block<<"*____STATISTICAL_PARAMETER_____ _______VALUE(S)_______"<<std::endl;
      *block<<" NUMBER OF OBSERVATIONS         "<<info.observationCount%"%22i"s<<std::endl;
      *block<<" NUMBER OF UNKNOWNS             "<<countParameter%"%22i"s<<std::endl;
      *block<<" NUMBER OF DEGREES OF FREEDOM   "<<(info.observationCount-countParameter)%"%22i"s<<std::endl;
      *block<<" WEIGHTED SQUARE SUM OF O-C     "<<info.lPl(0)%"%22.15e"s<<std::endl;
      if(sampling > 0)
        *block<<" SAMPLING INTERVAL (SECONDS)    "<<sampling%"%22f"s<<std::endl;
    }

    // ==================================================
    auto writeVector = [&](SinexBlockPtr block, const Vector x, const Vector sigma=Vector())
    {
      for(UInt i=0; i<x.size(); i++)
      {
        const std::string &object   = info.parameterName.at(i).object;
        const std::string &type     = info.parameterName.at(i).type;
        const std::string &temporal = info.parameterName.at(i).temporal;

        std::string parameterType;
        std::string siteCode       = "----";
        std::string pointCode      = "--";
        std::string unit           = "    ";
        std::string constraintCode = parameterIsConstrained.at(i) ? "1" : "2";

        // STAX, STAY, STAZ
        const Bool isStation = std::find(stationList.begin(), stationList.end(), object) != stationList.end();
        if(isStation && String::startsWith(type, "position."))
        {
          parameterType = "STA"+String::upperCase(type.substr(9, 1))+"  ";
          unit          = "m   ";
          siteCode      = String::upperCase(object);
          pointCode     = " A";
        }
        // XPO, YPO
        else if(object == "earth" && String::startsWith(type, "polarMotion.") && temporal.empty())
        {
          parameterType = String::upperCase(type.substr(12, 1)) + "PO   ";
          unit          = "mas ";
        }
        // XPOR, YPOR
        else if(object == "earth" && String::startsWith(type, "polarMotion.") && String::startsWith(temporal, "trend."))
        {
          parameterType = String::upperCase(type.substr(12, 1)) + "POR  ";
          unit          = "ma/d";
        }
        // UT1
        else if(object == "earth" && type == "UT1" && temporal.empty())
        {
          parameterType  = "UT    ";
          unit           = "ms  ";
          constraintCode = parameterIsConstrained.at(i) ? "0" : "2";
        }
        // LOD
        else if(object == "earth" && type == "UT1" && String::startsWith(temporal, "trend."))
        {
          parameterType = "LOD   ";
          unit          = "ms  ";
        }
        // SATA_X, SATA_Y, SATA_Z
        else if(!isStation && type.size() == 22 && type.substr(0, 14) == "antennaCenter.")
        {
          const std::string frequency = type.substr(type.size()-6, 2);
          if(frequency.at(0) != 'L')
            throw(Exception("unsupported antenna center parameter: " + info.parameterName.at(i).str()));
          if(type.at(14) == 'x')      parameterType  = "SATA_Y"; // swap X and Y names (definition different for GROOPS and IGS)
          else if(type.at(14) == 'y') parameterType  = "SATA_X"; // swap X and Y names (definition different for GROOPS and IGS)
          else if(type.at(14) == 'z') parameterType  = "SATA_Z";
          else
            throw(Exception("unsupported antenna center parameter: "+info.parameterName.at(i).str()));
          unit           = "m   ";
          siteCode       = object.substr(object.find('|')+1, 4); // SVN
          pointCode      = (frequency == "L*" ? "LC" : frequency);
          constraintCode = parameterIsConstrained.at(i) ? "0" : "2";
        }
        else
          throw(Exception("unknown parameter type: " + info.parameterName.at(i).str()));

        *block<<(i+1)%" %5i "s<<parameterType<<" "<<siteCode<<" "<<pointCode<<"    1 "<<Sinex::time2str(timeRef)
              <<" "<<unit<<" "<<constraintCode<<constraintCode<<x(i)%" %21.14e"s;
        if(sigma.size())
          *block<<sigma(i)%" %11.4e"s;
        *block<<std::endl;
      }
    };
    // ==================================================

    // add apriori antennaOffset of SATA_X, SATA_Y, SATA_Z
    // ---------------------------------------------------
    Vector xAprioriAntennaOffset(x.size());
    for(UInt i=0; i<info.parameterName.size(); i++)
    {
      const ParameterName &param = info.parameterName.at(i);
      if((param.type.size() == 22) && (param.type.substr(0, 14) == "antennaCenter."))
      {
        std::shared_ptr<PlatformGnssAntenna> antInfo;
        for(const auto &transmitterInfo : transmitterInfos)
        {
          antInfo = transmitterInfo.findEquipment<PlatformGnssAntenna>(timeRef);
          if(antInfo && (antInfo->str() == param.object))
            break;
        }
        if(!antInfo)
          throw(Exception("no antenna found for "+param.str()));
        const UInt idPattern = antInfo->antennaDef->findAntennaPattern(GnssType(param.type.substr(16, 6)), GnssAntennaDefinition::NoPatternFoundAction::USE_NEAREST_FREQUENCY);
        if(idPattern == NULLINDEX)
          throw(Exception("no antenna pattern found for "+param.str()));
        const Vector3d offset = antInfo->antennaDef->patterns.at(idPattern).offset + antInfo->local2antennaFrame.transform(antInfo->position);
        if(param.type.at(14) == 'x') xAprioriAntennaOffset(i) += offset.x();
        if(param.type.at(14) == 'y') xAprioriAntennaOffset(i) += offset.y();
        if(param.type.at(14) == 'z') xAprioriAntennaOffset(i) += offset.z();
      }
    }

    // -----------------
    if(x.size())
    {
      SinexBlockPtr block = sinex.addBlock("SOLUTION/ESTIMATE");
      *block<<"*INDEX _TYPE_ CODE PT SOLN _REF_EPOCH__ UNIT S ___ESTIMATED_VALUE___ __STD_DEV__"<<std::endl;
      writeVector(block, x+xAprioriAntennaOffset, sigmax.size() ? sigmax : Vector(x.size()));
    }
    // -----------------
    if(x0.size())
    {
      SinexBlockPtr block = sinex.addBlock("SOLUTION/APRIORI");
      *block<<"*INDEX _TYPE_ CODE PT SOLN _REF_EPOCH__ UNIT S ____APRIORI_VALUE____ __STD_DEV__"<<std::endl;
      writeVector(block, x0+xAprioriAntennaOffset, Vector(x0.size()));
    }
    // -----------------
    if(n.size())
    {
      SinexBlockPtr block = sinex.addBlock("SOLUTION/NORMAL_EQUATION_VECTOR");
      *block<<"*INDEX _TYPE_ CODE PT SOLN _REF_EPOCH__ UNIT S ___RIGHT_HAND_SIDE___"<<std::endl;
      writeVector(block, n);
    }

    // ==================================================
    auto writeMatrix = [&](SinexBlockPtr block, const Matrix &N)
    {
      for(UInt i=0; i<N.rows(); i++)
        for(UInt k=i; k<N.rows(); k++)
          if(N(i,k))
          {
            *block<<(i+1)%" %5i"s<<(k+1)%" %5i"s<<N(i, k)%" %21.14e"s;
            for(UInt l=1; (l<3) && (k+1<N.rows()) && N(i,k+1); l++, k++)
              *block<<N(i, k+1)%" %21.14e"s;
            *block<<std::endl;
          }
    };
    // ==================================================

    if(N.size())
    {
      SinexBlockPtr block = sinex.addBlock("SOLUTION/NORMAL_EQUATION_MATRIX U");
      *block<<"*PARA1 PARA2 _______PARA2+0_______ _______PARA2+1_______ _______PARA2+2_______"<<std::endl;
      writeMatrix(block, N);
    }
    // -----------------
    if(dN.size())
    {
      SinexBlockPtr block = sinex.addBlock("SOLUTION/MATRIX_APRIORI U INFO");
      *block<<"*PARA1 PARA2 _______PARA2+0_______ _______PARA2+1_______ _______PARA2+2_______"<<std::endl;
      writeMatrix(block, dN);
    }

    // ==================================================

    // write SINEX files
    // -----------------
    if(!fileNameSinexNormals.empty())
    {
      logStatus<<"write full SINEX file <"<<fileNameSinexNormals<<">"<<Log::endl;
      sinex.header.replace(60, 5, countParameter%"%05i"s);
      writeFileSinex(fileNameSinexNormals, sinex);
    }
    if(!fileNameSinexCoords.empty())
    {
      logStatus<<"write coordinates SINEX file <"<<fileNameSinexCoords<<">"<<Log::endl;
      sinex.blocks.remove_if([](const SinexBlockPtr &b){return b->label == "SOLUTION/NORMAL_EQUATION_VECTOR";});
      sinex.blocks.remove_if([](const SinexBlockPtr &b){return b->label == "SOLUTION/NORMAL_EQUATION_MATRIX U";});
      sinex.blocks.remove_if([](const SinexBlockPtr &b){return b->label == "SOLUTION/MATRIX_APRIORI U INFO";});
      sinex.header.replace(60, 5, countParameter%"%05i"s);
      writeFileSinex(fileNameSinexCoords, sinex);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
