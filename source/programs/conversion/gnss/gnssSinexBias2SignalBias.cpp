/***********************************************/
/**
* @file gnssSinexBias2SignalBias.cpp
*
* @brief Converts GNSS signal biases from SINEX Bias format to GnssSignalBias format.
*
* @author Sebastian Strasser
* @date 2020-09-30
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Converts GNSS signal biases from \href{https://files.igs.org/pub/data/format/sinex_bias_100.pdf}{IGS SINEX Bias format}
to \file{GnssSignalBias format}{gnssSignalBias}.

Only satellite observable-specific signal biases (OSB) are supported at the moment.
If multiple entries exist for the same bias, the weighted average (based on time span) of all entries is used.
Time-variable biases are not supported at the moment.

See also \program{GnssSignalBias2SinexBias}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "files/fileGnssSignalBias.h"
#include "files/fileGnssReceiverDefinition.h"
#include "inputOutput/fileSinex.h"

/***** CLASS ***********************************/

/** @brief Converts GNSS signal biases from SINEX Bias format to GnssSignalBias format.
* @ingroup programsGroup */
class GnssSinexBias2SignalBias
{
  class Bias
  {
  public:
    Time timeStart, timeEnd;
    Double value, slope;
  };

public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GnssSinexBias2SignalBias, SINGLEPROCESS, "Converts GNSS signal biases from SINEX Bias format to GnssSignalBias format.", Conversion, Gnss)

/***********************************************/

void GnssSinexBias2SignalBias::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameOutSignalBias, fileNameInSinexBias, fileNameInGlonassReceiverDefinition;
    std::vector<std::string> identifiers;

    readConfig(config, "outputfileSignalBias",               fileNameOutSignalBias,               Config::MUSTSET,  "", "identifier is appended to file name");
    readConfig(config, "inputfileSinexBias",                 fileNameInSinexBias,                 Config::MUSTSET,  "", "");
    readConfig(config, "inputfileGlonassReceiverDefinition", fileNameInGlonassReceiverDefinition, Config::OPTIONAL, "{groopsDataDir}/gnss/transmitterGlonass/receiverDefinition/receiverDefinition.xml", "GLONASS frequency number mapping");
    readConfig(config, "identifier",                         identifiers,                         Config::OPTIONAL, "", "(empty = all) satellite PRN, e.g. G23 or E05");
    if(isCreateSchema(config)) return;

    logStatus<<"read SINEX Bias file <"<<fileNameInSinexBias<<">"<<Log::endl;
    Sinex sinexFile(fileNameInSinexBias);
    std::vector<std::string> lines = sinexFile.getBlock<Sinex::SinexText>("BIAS/SOLUTION")->lines();

    std::vector<GnssReceiverDefinitionPtr> receiverDefinitions;
    if(!fileNameInGlonassReceiverDefinition.empty())
    {
      logStatus<<"read GLONASS receiver definition file <"<<fileNameInGlonassReceiverDefinition<<">"<<Log::endl;
      readFileGnssReceiverDefinition(fileNameInGlonassReceiverDefinition, receiverDefinitions);
    }

    Bool printStationWarning = TRUE;
    Bool printNonOsbWarning = TRUE;
    std::map<std::string, std::map<GnssType, std::vector<Bias>>> id2type2biases;
    for(const auto &line : lines)
    {
      std::string biasType = String::trim(line.substr(1,3));
      std::string svn      = String::trim(line.substr(6,4));
      std::string prn      = String::trim(line.substr(11,3));
      std::string station  = String::trim(line.substr(15,9));

      if(biasType != "OSB")
      {
        if(printNonOsbWarning)
        {
          logWarning<<"bias types other than OSB (observable-specific signal bias) are not yet implemented and will be ignored."<<Log::endl;
          printNonOsbWarning = FALSE;
        }
        continue;
      }
      if(station.size() && (!identifiers.size() || std::find(identifiers.begin(), identifiers.end(), station) != identifiers.end()))
      {
        if(printStationWarning)
        {
          logWarning<<"station biases are not yet implemented and will be ignored."<<Log::endl;
          printStationWarning = FALSE;
        }
        continue;
      }
      if(identifiers.size() && std::find(identifiers.begin(), identifiers.end(), prn) == identifiers.end())
        continue;

      GnssType gnssType = GnssType(line.substr(25,3)+prn);
      if(gnssType == GnssType::GLONASS)
      {
        auto iter = std::find_if(receiverDefinitions.begin(), receiverDefinitions.end(), [&](const GnssReceiverDefinitionPtr &def){ return def->serial == svn; });
        if(iter != receiverDefinitions.end())
          gnssType.setFrequencyNumber(String::toInt((*iter)->version));
      }
      std::string unit = String::trim(line.substr(65,4));
      Bias bias;
      bias.value = String::toDouble(line.substr(70,21));
      bias.slope = line.size() >= 125 ? String::toDouble(line.substr(104,21)) : 0.;
      bias.timeStart = Sinex::str2time(line, 35, TRUE);
      bias.timeEnd   = Sinex::str2time(line, 50, TRUE);
      if(unit == "ns")
      {
        bias.value *= 1e-9*LIGHT_VELOCITY;
        bias.slope *= 1e-9*LIGHT_VELOCITY;
      }
      else if(unit == "cyc")
      {
        bias.value *= gnssType.wavelength();
        bias.slope *= gnssType.wavelength();
      }
      else
        throw(Exception("unknown unit: "+unit));
      id2type2biases[prn][gnssType].push_back(bias);
    }

    if(id2type2biases.size())
    {
      logStatus<<"write signal bias files <"<<fileNameOutSignalBias<<">"<<Log::endl;
      for(const auto &id : id2type2biases)
      {
        GnssSignalBias signalBias;
        for(const auto &type : id.second)
        {
          signalBias.type.push_back(type.first);
          Double weightedSum  = std::accumulate(type.second.begin(), type.second.end(), 0.,
                                                [](Double a, const Bias &bias){ return a + bias.value*(bias.timeEnd-bias.timeStart).seconds(); });
          Double sumOfWeights = std::accumulate(type.second.begin(), type.second.end(), 0.,
                                                [](Double a, const Bias &bias){ return a + (bias.timeEnd-bias.timeStart).seconds(); });
          signalBias.bias.push_back(weightedSum/sumOfWeights);
        }
        writeFileGnssSignalBias(fileNameOutSignalBias.appendBaseName("."+id.first), signalBias);
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
