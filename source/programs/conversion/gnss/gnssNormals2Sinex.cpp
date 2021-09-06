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
Metadata input files include \configFile{stationInfo/transmitterInfo}{gnssStationInfo}, \configFile{antennaDefinition}{gnssAntennaDefinition},
and \configFile{stationList/transmitterList}{stringList}, see \program{GnssAntex2AntennaDefinition}.

See also \program{Sinex2Normals} and \program{NormalsSphericalHarmonics2Sinex}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "inputOutput/fileSinex.h"
#include "files/fileGnssStationInfo.h"
#include "files/fileMatrix.h"
#include "files/fileNormalEquation.h"
#include "files/fileStringTable.h"

/***** CLASS ***********************************/

/** @brief Write GNSS data/metadata and normal equations to SINEX format.
* @ingroup programsConversionGroup */
class GnssNormals2Sinex
{
  static void addVector(Sinex::SinexSolutionVectorPtr vector, const Time &time, const std::vector<ParameterName> &parameterName, const Vector x,
                        const Vector sigma, const std::vector<Bool> &parameterIsContrained, const std::vector<std::string> &stationList);
  static void addAprioriAntennaOffset(const std::map<std::string, GnssAntennaInfo> &antennas, Bool addEccentricity, Bool swapXY, const std::vector<ParameterName> &parameterNames, MatrixSliceRef x);

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
  readConfig(config, "inputfileTransmitterInfo",   var.fileNameTransmitterInfo, Config::MUSTSET,  "{groopsDataDir}/gnss/transmitter/transmitterInfo/igs/igs14/transmitterInfo_igs14.{prn}.xml", "transmitter info file template");
  readConfig(config, "inputfileAntennaDefinition", var.fileNameAntennaDef,      Config::MUSTSET,  "{groopsDataDir}/gnss/transmitter/antennaDefinition/igs/igs14/transmitterInfo_igs14.dat",     "transmitter phase centers and variations (ANTEX)");
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
    Sinex       sinex;
    std::string variableStationName, antennaModel;
    Double      sampling = 0;
    std::vector<TransmitterConstellation> constellations;

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
      readConfig(config, "inputfileAntennaDefinition", fileNameAntennaDef,  Config::MUSTSET,  "{groopsDataDir}/gnss/receiverStation/antennaDefinition/igs/igs14/antennaDefinition_igs14.dat", "station phase centers and variations (ANTEX)");
      readConfig(config, "variableStationName",        variableStationName, Config::DEFAULT,  "station", "loop variable for station names from station list");
      readConfig(config, "observationTimeStart",       timeStartObs,        Config::MUSTSET,  "", "start time for which solution has observations");
      readConfig(config, "observationTimeEnd",         timeEndObs,          Config::MUSTSET,  "", "end time for which solution has observations");
      endSequence(config);
    }
    readConfig(config, "time",                       timeRef,             Config::MUSTSET,  "", "reference time for parameters");
    readConfig(config, "sampling",                   sampling,            Config::OPTIONAL, "", "[seconds] observation sampling");
    readConfig(config, "antennaCalibrationModel",    antennaModel,        Config::MUSTSET,  "", "e.g. IGS14_WWWW (WWWW = ANTEX release GPS week)");
    sinex.readConfigHeader(config);
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
    const UInt countParameter = N.rows();

    Matrix dN;
    std::vector<Bool> parameterIsConstrained(countParameter, FALSE);
    if(!fileNameAprMat.empty())
    {
      logStatus<<"reading normal equation matrix of applied constraints <"<<fileNameAprMat<<">"<<Log::endl;
      Vector n;
      NormalEquationInfo info;
      readFileNormalEquation(fileNameAprMat, info, dN, n);
      if(dN.rows() != parameterIsConstrained.size())
        throw(Exception("Parameter count in constraint matrix and normal equation matrix differs (" + dN.rows()%"%i"s + " vs. "+ N.rows()%"%i"s +" )."));
      for(UInt i = 0; i < dN.rows(); i++)
        parameterIsConstrained.at(i) = (dN(i, i) != 0.0);
    }

    logStatus<<"reading station list from <"<<fileNameStationList<<">"<<Log::endl;
    std::vector<std::string> stationList;
    readFileStringList(fileNameStationList, stationList);

    logStatus<<"reading station antenna definitions from <"<<fileNameAntennaDef<<">"<<Log::endl;
    std::vector<GnssAntennaDefinitionPtr> antennaDefinitionList;
    readFileGnssAntennaDefinition(fileNameAntennaDef, antennaDefinitionList);

    logStatus<<"reading station infos from <"<<fileNameStationInfo<<">"<<Log::endl;
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

    std::vector<GnssStationInfo> transmitterInfos;
    for(const auto &constellation : constellations)
    {
      logStatus<<"reading transmitter list from <"<<constellation.fileNameTransmitterList<<">"<<Log::endl;
      std::vector<std::string> transmitterList;
      readFileStringList(constellation.fileNameTransmitterList, transmitterList);

      logStatus<<"reading transmitter antenna definitions from <"<<constellation.fileNameAntennaDef<<">"<<Log::endl;
      std::vector<GnssAntennaDefinitionPtr> antennaDefinitionList;
      readFileGnssAntennaDefinition(constellation.fileNameAntennaDef, antennaDefinitionList);

      logStatus<<"reading transmitter infos from <"<<constellation.fileNameTransmitterInfo<<">"<<Log::endl;
      addVariable(constellation.variablePrn, fileNameVariableList);
      for(const auto &prn : transmitterList)
      {
        fileNameVariableList[constellation.variablePrn]->setValue(prn);
        GnssStationInfo transmitterInfo;
        readFileGnssStationInfo(constellation.fileNameTransmitterInfo(fileNameVariableList), transmitterInfo);
        transmitterInfo.fillAntennaPattern(antennaDefinitionList);
        transmitterInfos.push_back(transmitterInfo);
      }
    }

    // ==================================================

    // add data to SINEX
    // -----------------
    sinex.addSiteIdBlock(stationInfos);
    sinex.addSiteReceiverBlock(stationInfos, timeRef);
    sinex.addSiteAntennaBlock(stationInfos, timeRef);
    sinex.addSiteGpsPhaseCenterBlock(stationInfos, timeRef, antennaModel);
    sinex.addSiteGalileoPhaseCenterBlock(stationInfos, timeRef, antennaModel);
    sinex.addSiteEccentricityBlock(stationInfos, timeRef);
    sinex.addSatelliteIdBlock(transmitterInfos, timeRef);
    sinex.addSatellitePhaseCenter(transmitterInfos, timeRef, antennaModel);

    // SOLUTION/EPOCHS
    Sinex::SinexSolutionEpochsPtr solutionEpochs = sinex.addBlock<Sinex::SinexSolutionEpochs>("SOLUTION/EPOCHS");
    for(const auto &stationInfo : stationInfos)
      solutionEpochs->addEpoch(Sinex::Epoch(String::upperCase(stationInfo.markerName), "A", "1", "P", timeStartObs, timeEndObs));

    // SOLUTION/STATISTICS
    Sinex::SinexSolutionStatisticsPtr solutionStatistics = sinex.addBlock<Sinex::SinexSolutionStatistics>("SOLUTION/STATISTICS");
    solutionStatistics->addValue("NUMBER OF OBSERVATIONS", info.observationCount);
    solutionStatistics->addValue("NUMBER OF UNKNOWNS", countParameter);
    solutionStatistics->addValue("NUMBER OF DEGREES OF FREEDOM", info.observationCount-countParameter);
    solutionStatistics->addValue("WEIGHTED SQUARE SUM OF O-C", info.lPl(0));
    if(sampling > 0)
      solutionStatistics->addValue("SAMPLING INTERVAL (SECONDS)", sampling);

    auto getTransmitterName = [](const GnssStationInfo &info, UInt idAnt)     { return info.antenna.at(idAnt).name+"|"+info.antenna.at(idAnt).serial+"|"+info.antenna.at(idAnt).radome; };
    auto getStationName     = [](const GnssStationInfo &info, UInt /*idAnt*/) { return String::lowerCase(info.markerName); };
    auto getAntennas = [&](const std::vector<GnssStationInfo> &infos, std::function<std::string(const GnssStationInfo&, UInt)> getName)
    {
      std::map<std::string, GnssAntennaInfo> antennas;
      for(const auto &info : infos)
      {
        const UInt idAnt = info.findAntenna(timeRef);
        if(idAnt == NULLINDEX)
        {
          logWarning << info.markerName << ": no antenna found at " << timeRef.dateTimeStr() << Log::endl;
          continue;
        }
        antennas[getName(info, idAnt)] = info.antenna.at(idAnt);
      }
      return antennas;
    };

    // SOLUTION/ESTIMATE
    if(x.size())
    {
      addAprioriAntennaOffset(getAntennas(transmitterInfos, getTransmitterName), TRUE/*addEccentricity*/,  TRUE/*swapXY*/,  info.parameterName, x);
      addAprioriAntennaOffset(getAntennas(stationInfos,     getStationName),     FALSE/*addEccentricity*/, FALSE/*swapXY*/, info.parameterName, x);

      Sinex::SinexSolutionVectorPtr solutionEstimate = sinex.addBlock<Sinex::SinexSolutionVector>("SOLUTION/ESTIMATE");
      addVector(solutionEstimate, timeRef, info.parameterName, x, sigmax.size() ? sigmax : Vector(), parameterIsConstrained, stationList);
    }

    // SOLUTION/APRIORI
    if(x0.size())
    {
      addAprioriAntennaOffset(getAntennas(transmitterInfos, getTransmitterName), TRUE/*addEccentricity*/,  TRUE/*swapXY*/,  info.parameterName, x0);
      addAprioriAntennaOffset(getAntennas(stationInfos,     getStationName),     FALSE/*addEccentricity*/, FALSE/*swapXY*/, info.parameterName, x0);

      for(UInt i = 0; i < x0.size(); i++)
        if(dN.size() && !parameterIsConstrained.at(i))
          sigmax0(i) = 0;
      Sinex::SinexSolutionVectorPtr solutionApriori = sinex.addBlock<Sinex::SinexSolutionVector>("SOLUTION/APRIORI");
      addVector(solutionApriori, timeRef, info.parameterName, x0, sigmax0, parameterIsConstrained, stationList);
    }

    // SOLUTION/NORMAL_EQUATION_VECTOR
    Sinex::SinexSolutionVectorPtr solutionNormalEquationVector = sinex.addBlock<Sinex::SinexSolutionVector>("SOLUTION/NORMAL_EQUATION_VECTOR");
    addVector(solutionNormalEquationVector, timeRef, info.parameterName, n, Vector(), parameterIsConstrained, stationList);

    // SOLUTION/NORMAL_EQUATION_MATRIX
    Sinex::SinexSolutionMatrixPtr solutionNormalEquationMatrix = sinex.addBlock<Sinex::SinexSolutionMatrix>("SOLUTION/NORMAL_EQUATION_MATRIX "s + (N.isUpper() ? "U" : "L"));
    solutionNormalEquationMatrix->setMatrix(N);

    // SOLUTION/MATRIX_APRIORI
    if(dN.size())
    {
      Sinex::SinexSolutionMatrixPtr solutionNormalAprioriMatrix = sinex.addBlock<Sinex::SinexSolutionMatrix>("SOLUTION/MATRIX_APRIORI "s + (dN.isUpper() ? "U" : "L") + " INFO");
      solutionNormalAprioriMatrix->setMatrix(dN);
    }

    // ==================================================

    // write SINEX files
    // -----------------
    if(!fileNameSinexNormals.empty())
    {
      logStatus<<"write full SINEX file <"<<fileNameSinexNormals<<">"<<Log::endl;
      sinex.writeFile(fileNameSinexNormals);
    }
    if(!fileNameSinexCoords.empty())
    {
      logStatus<<"write coordinates SINEX file <"<<fileNameSinexCoords<<">"<<Log::endl;
      sinex.removeBlock("SOLUTION/NORMAL_EQUATION_VECTOR");
      sinex.removeBlock("SOLUTION/NORMAL_EQUATION_MATRIX "s + (N.isUpper() ? "U" : "L"));
      if(dN.size())
        sinex.removeBlock("SOLUTION/MATRIX_APRIORI "s + (dN.isUpper() ? "U" : "L") + " INFO");
      sinex.writeFile(fileNameSinexCoords);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssNormals2Sinex::addVector(Sinex::SinexSolutionVectorPtr vector, const Time &time0, const std::vector<ParameterName> &parameterName, const Vector x,
                                  const Vector sigma, const std::vector<Bool> &parameterIsContrained, const std::vector<std::string> &stationList)
{
  try
  {
    for(UInt i = 0; i < x.size(); i++)
    {
      const std::string object   = parameterName.at(i).object;
      const std::string type     = parameterName.at(i).type;
      const std::string temporal = parameterName.at(i).temporal;

      Sinex::Parameter parameter;
      parameter.parameterIndex = i+1;
      parameter.solutionId     = "1";
      parameter.constraintCode = parameterIsContrained.at(i) ? "1" : "2";
      parameter.time           = time0;
      parameter.value          = x(i);
      if(sigma.size())
        parameter.sigma        = sigma(i);

      // STAX, STAY, STAZ
      const Bool isStation = std::find(stationList.begin(), stationList.end(), object) != stationList.end();
      if(isStation && String::startsWith(type, "position."))
      {
        parameter.parameterType = "STA" + String::upperCase(type.substr(9, 1));
        parameter.unit          = "m";
        parameter.siteCode      = String::upperCase(object);
        parameter.pointCode     = "A";
      }
      // XPO, YPO
      else if(object == "earth" && String::startsWith(type, "polarMotion.") && temporal.empty())
      {
        parameter.parameterType = String::upperCase(type.substr(12, 1)) + "PO";
        parameter.unit          = "mas";
      }
      // XPOR, YPOR
      else if(object == "earth" && String::startsWith(type, "polarMotion.") && String::startsWith(temporal, "trend."))
      {
        parameter.parameterType = String::upperCase(type.substr(12, 1)) + "POR";
        parameter.unit          = "ma/d";
      }
      // UT1
      else if(object == "earth" && type == "UT1" && temporal.empty())
      {
        parameter.parameterType  = "UT";
        parameter.unit           = "ms";
        parameter.constraintCode = parameterIsContrained.at(i) ? "0" : "2";
      }
      // LOD
      else if(object == "earth" && type == "UT1" && String::startsWith(temporal, "trend."))
      {
        parameter.parameterType = "LOD";
        parameter.unit          = "ms";
      }
      // SATA_X, SATA_Y, SATA_Z
      else if(!isStation && type.size() == 22 && type.substr(0, 14) == "antennaCenter.")
      {
        const std::string frequency = type.substr(type.size()-6, 2);
        if(frequency.at(0) != 'L')
          throw(Exception("unsupported antenna center parameter: " + parameterName.at(i).str()));

        if(     type.at(14) == 'x') parameter.parameterType  = "SATA_Y"; // swap X and Y names (definition different for GROOPS and IGS)
        else if(type.at(14) == 'y') parameter.parameterType  = "SATA_X"; // swap X and Y names (definition different for GROOPS and IGS)
        else if(type.at(14) == 'z') parameter.parameterType  = "SATA_Z";
        else throw(Exception("unsupported antenna center parameter: " + parameterName.at(i).str()));
        parameter.unit           = "m";
        parameter.siteCode       = object.substr(object.find('|')+1, 4); // SVN
        parameter.pointCode      = (frequency == "L*" ? "LC" : frequency);
        parameter.constraintCode = parameterIsContrained.at(i) ? "0" : "2";
      }
      else
        throw(Exception("unknown parameter type: " + parameterName.at(i).str()));

      vector->addParameter(parameter);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssNormals2Sinex::addAprioriAntennaOffset(const std::map<std::string, GnssAntennaInfo> &antennas, Bool addEccentricity, Bool swapXY, const std::vector<ParameterName> &parameterNames, MatrixSliceRef x)
{
  try
  {
    for(UInt i = 0; i < x.size(); i++)
    {
      const ParameterName* param = &parameterNames.at(i);

      if(param->type.size() == 22 && param->type.substr(0, 14) == "antennaCenter." && antennas.find(param->object) != antennas.end())
      {
        const GnssType type(param->type.substr(16, 6));
        const UInt idPattern = antennas.at(param->object).antennaDef->findAntennaPattern(type, GnssAntennaDefinition::NoPatternFoundAction::USE_NEAREST_FREQUENCY);
        if(idPattern == NULLINDEX)
          throw(Exception("no antenna pattern found for " + param->str()));

        const Vector3d offset = antennas.at(param->object).local2antennaFrame.inverseTransform(antennas.at(param->object).antennaDef->pattern.at(idPattern).offset);
        if(param->type.at(14) == (swapXY ? 'y' : 'x'))
          x(i,0) += (addEccentricity ? offset.x() + antennas.at(param->object).position.x() : offset.x());
        if(param->type.at(14) == (swapXY ? 'x' : 'y'))
          x(i,0) += (addEccentricity ? offset.y() + antennas.at(param->object).position.y() : offset.y());
        if(param->type.at(14) == 'z')
          x(i,0) += (addEccentricity ? offset.z() + antennas.at(param->object).position.z() : offset.z());
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
