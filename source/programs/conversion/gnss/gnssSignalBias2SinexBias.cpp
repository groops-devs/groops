/***********************************************/
/**
* @file gnssSignalBias2SinexBias.cpp
*
* @brief Convert GNSS signal biases from GROOPS format to IGS SINEX Bias format.
**
* @author Sebastian Strasser
* @date 2019-01-24
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Convert \file{GNSS signal biases}{gnssSignalBias} from GROOPS format to \href{https://files.igs.org/pub/data/format/sinex_bias_100.pdf}{IGS SINEX Bias format}.
Biases can be provided via \config{transmitterBiases} and/or \config{receiverBiases}.
Phase biases without attribute (e.g. \verb|L1*|) are automatically expanded so each code
bias has a corresponding phase bias
(Example: \verb|C1C|, \verb|C1W|, \verb|L1*| are converted to \verb|C1C|, \verb|C1W|, \verb|L1C|, \verb|L1W|).

Time-variable biases (e.g. GPS L5 satellite phase bias) can be provided via \config{timeVariableBias}.
Their time span will be based on the provided epochs ($t \pm \Delta t / 2$).
The slope of the bias can be optionally provided in the second data column.

If GLONASS receiver biases depend on frequency number, those must be defined in \configFile{inputfileTransmitterInfo}{gnssStationInfo}
to get the correct PRN/SVN assignment to the biases.

See IGS SINEX Bias format description for further details on header information.

See also \program{GnssSinexBias2SignalBias} and \program{GnssBiasClockAlignment}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileGnssSignalBias.h"
#include "files/fileGnssStationInfo.h"
#include "files/fileMatrix.h"
#include "inputOutput/file.h"
#include "inputOutput/fileSinex.h"
#include "inputOutput/system.h"

/***** CLASS ***********************************/

/** @brief Convert GNSS signal biases from GROOPS format to IGS SINEX Bias format.
* @ingroup programsConversionGroup */
class GnssSignalBias2SinexBias
{
public:
  class TimeVariableBias
  {
  public:
    FileName inNameBias;
    GnssType type;
    Matrix   biases;
  };

  class Data
  {
  public:
    FileName inNameBias, inNameTransmitterInfo;
    std::string identifier;
    GnssSignalBias biases;
    std::vector<TimeVariableBias> timeVariableBiases;

    inline Bool isTimeVariableBias(const GnssType &type) const
    {
      return std::find_if(timeVariableBiases.begin(), timeVariableBiases.end(), [&](const TimeVariableBias &bias){ return bias.type.type == (type & GnssType::NOPRN).type; }) != timeVariableBiases.end();
    }
  };

private:
  std::map<std::string, GnssStationInfo> prn2info;
  std::map<Int, std::vector<std::string>> freqNo2prns;

  void readData(std::vector<Data> &data) const;
  void writeLine(OutFile &file, const Time &timeStart, const Time &timeEnd, std::string prn, std::string svn, std::string stationName,
                 const GnssType &type, std::string unit, Double bias, Double sigma, Double biasSlope, Double biasSlopeSigma) const;
  void writeData(OutFile &file, const Time &timeStart, const Time &timeEnd, const Data &data, Bool isStation) const;
  UInt countBiases(const std::vector<Data> &data) const;

public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(GnssSignalBias2SinexBias, SINGLEPROCESS, "Convert GNSS signal biases from GROOPS format to IGS SINEX Bias format.", Conversion, Gnss)

/***********************************************/

template<> Bool readConfig(Config &config, const std::string &name, GnssSignalBias2SinexBias::TimeVariableBias &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  if(!readConfigSequence(config, name, mustSet, defaultValue, annotation))
    return FALSE;
  readConfig(config, "inputfileSignalBias",  var.inNameBias,  Config::MUSTSET, "", "columns: mjd, bias [m], (biasSlope [m/s])");
  readConfig(config, "type",                 var.type,        Config::MUSTSET, "", "bias type");
  endSequence(config);
  return TRUE;
}

/***********************************************/

template<> Bool readConfig(Config &config, const std::string &name, GnssSignalBias2SinexBias::Data &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  if(!readConfigSequence(config, name, mustSet, defaultValue, annotation))
    return FALSE;
  readConfig(config, "inputfileSignalBias",  var.inNameBias,         Config::MUSTSET,  "", "signal bias file");
  readConfig(config, "timeVariableBias",     var.timeVariableBiases, Config::OPTIONAL, "", "one entry per time variable bias type");
  readConfig(config, "identifier",           var.identifier,         Config::MUSTSET,  "", "PRN or station name (e.g. G23 or wtzz)");
  endSequence(config);
  return TRUE;
}

/***********************************************/

void GnssSignalBias2SinexBias::readData(std::vector<Data> &data) const
{
  try
  {
    auto iter = data.begin();
    while(iter != data.end())
    {
      try
      {
        readFileGnssSignalBias(iter->inNameBias, iter->biases);
        for(auto &type : iter->biases.type)
          if(type == (GnssType::GALILEO + GnssType::UNKNOWN_ATTRIBUTE))
            type = (type & ~GnssType::ATTRIBUTE) + GnssType::X;
        for(auto &&bias : iter->timeVariableBiases)
        {
          try
          {
            readFileMatrix(bias.inNameBias, bias.biases);
          }
          catch(std::exception &e)
          {
            logWarning << e.what() << " continue..." << Log::endl;
            continue;
          }
        }
      }
      catch(std::exception &e)
      {
        logWarning << e.what() << " continue..." << Log::endl;
        iter = data.erase(iter);
        continue;
      }

      iter++;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssSignalBias2SinexBias::writeLine(OutFile &file, const Time &timeStart, const Time &timeEnd, std::string prn, std::string svn, std::string stationName,
                                         const GnssType &type, std::string unit, Double bias, Double biasSigma, Double biasSlope, Double biasSlopeSigma) const
{
  std::transform(stationName.begin(), stationName.end(), stationName.begin(), ::toupper);
  prn.resize(3, ' ');
  svn.resize(4, ' ');
  unit.resize(4, ' ');
  stationName.resize(9, ' ');

  if(stationName.find_first_not_of(' ') != std::string::npos)
  {
    if(prn.find_first_not_of(' ') == std::string::npos)
      prn.at(0) = type.str().at(3);
    if(svn.find_first_not_of(' ') == std::string::npos)
      svn.at(0) = type.str().at(3);
  }

  file << " OSB  " << svn << ' ' << prn << ' ' << stationName << ' ' << type.str().substr(0,3) << std::string(7, ' ')
       << Sinex::time2str(timeStart, TRUE) << ' ' << Sinex::time2str(timeEnd, TRUE) << ' ' << unit << bias%" % 21.14e"s << biasSigma%" % 11.5e"s;
  if(biasSlope != 0.)
    file << biasSlope%" % 21.14e"s << biasSlopeSigma%" % 11.5e"s;
  file << std::endl;
}

/***********************************************/

void GnssSignalBias2SinexBias::writeData(OutFile &file, const Time &timeStart, const Time &timeEnd, const Data &data, Bool isStation) const
{
  try
  {
    auto getSvn = [&](const std::string &prn)
    {
      UInt idAnt = prn2info.at(prn).findAntenna(timeStart);
      if(idAnt == NULLINDEX)
        throw(Exception("no SVN found for " + data.identifier));

      return prn2info.at(prn).antenna.at(idAnt).serial;
    };

    std::string prn = (isStation ? "" : data.identifier);
    std::string svn = (isStation ? "" : getSvn(prn));
    std::string stationName = (isStation ? data.identifier : "");

    auto writeCodeBias = [&](const Time &timeStart, const Time &timeEnd, const GnssType &type, Double bias, Double biasSlope=0.)
    {
      // GLONASS station bias per satellite
      if(isStation && type == GnssType::GLONASS && type.frequencyNumber() != 9999)
      {
        for(const auto &prn : freqNo2prns.at(type.frequencyNumber()))
          writeLine(file, timeStart, timeEnd, prn, getSvn(prn), stationName, type, "ns", bias/LIGHT_VELOCITY*1e9, 0, biasSlope/LIGHT_VELOCITY*1e9, 0);
        return;
      }

      writeLine(file, timeStart, timeEnd, prn, svn, stationName, type, "ns", bias/LIGHT_VELOCITY*1e9, 0, biasSlope/LIGHT_VELOCITY*1e9, 0);
    };

    auto writePhaseBias = [&](const Time &timeStart, const Time &timeEnd, const GnssType &type, Double bias, Double biasSlope=0.)
    {
      for(const auto &type2 : data.biases.type)
        if(type2 == (type & ~GnssType::TYPE) + GnssType::RANGE)
        {
          // GLONASS station bias per satellite
          if(isStation && type == GnssType::GLONASS && type.frequencyNumber() != 9999)
          {
            for(const auto &prn : freqNo2prns.at(type.frequencyNumber()))
              writeLine(file, timeStart, timeEnd, prn, getSvn(prn), stationName, (type2 & ~GnssType::TYPE) + GnssType::PHASE,
                        "ns", bias/LIGHT_VELOCITY*1e9, 0, biasSlope/LIGHT_VELOCITY*1e9, 0);
            continue;
          }

          writeLine(file, timeStart, timeEnd, prn, svn, stationName, (type2 & ~GnssType::TYPE) + GnssType::PHASE,
                    "ns", bias/LIGHT_VELOCITY*1e9, 0, biasSlope/LIGHT_VELOCITY*1e9, 0);
        }
    };


    // Code biases
    for(UInt i = 0; i < data.biases.type.size(); i++)
      if(data.biases.type.at(i) == GnssType::RANGE && !data.isTimeVariableBias(data.biases.type.at(i)))
        writeCodeBias(timeStart, timeEnd, data.biases.type.at(i), data.biases.bias.at(i));

    // Phase biases (for wildcard types, i.e. L1*, write one line per matching code bias type with the same phase bias; example: C1C, C1W, L1* ==> C1C, C1W, L1C, L1W)
    for(UInt i = 0; i < data.biases.type.size(); i++)
      if(data.biases.type.at(i) == GnssType::PHASE && !data.isTimeVariableBias(data.biases.type.at(i)))
        writePhaseBias(timeStart, timeEnd, data.biases.type.at(i), data.biases.bias.at(i));

    // Time-variable biases
    for(const auto &bias : data.timeVariableBiases)
    {
      Double constBias = 0;
      for(UInt i = 0; i < data.biases.type.size(); i++)
        if((data.biases.type.at(i) & GnssType::NOPRN).type == bias.type.type) // exact match
          {
            constBias = data.biases.bias.at(i);
            break;
          }

      std::vector<Time> times(bias.biases.rows());
      for(UInt idEpoch = 0; idEpoch < bias.biases.rows(); idEpoch++)
        times.at(idEpoch) = mjd2time(bias.biases(idEpoch, 0));
      Time sampling = medianSampling(times);
      for(UInt idEpoch = 0; idEpoch < bias.biases.rows(); idEpoch++)
      {
        const Time timeStart = times.at(idEpoch)-0.5*sampling;
        const Time timeEnd   = times.at(idEpoch)+0.5*sampling;
        const Double biasValue = constBias + bias.biases(idEpoch, 1);
        const Double biasSlope = bias.biases.columns() > 2 ? bias.biases(idEpoch, 2) : 0;
        if(bias.type == GnssType::RANGE)
          writeCodeBias(timeStart, timeEnd, bias.type, biasValue, biasSlope);
        else if(bias.type == GnssType::PHASE)
          writePhaseBias(timeStart, timeEnd, bias.type, biasValue, biasSlope);
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

UInt GnssSignalBias2SinexBias::countBiases(const std::vector<Data> &data) const
{
  try
  {
    UInt count = 0;
    for(const auto &d : data)
    {
      for(UInt i = 0; i < d.biases.type.size(); i++)
        if(d.biases.type.at(i) == GnssType::RANGE)
          count++;

      // Phase biases (for wildcard types, i.e. L1*, write one line per matching code bias type with the same phase bias; example: C1C, C1W, L1* ==> C1C, C1W, L1C, L1W)
      for(UInt i = 0; i < d.biases.type.size(); i++)
        if(d.biases.type.at(i) == GnssType::PHASE)
          for(UInt j = 0; j < d.biases.type.size(); j++)
            if(d.biases.type.at(j) == (d.biases.type.at(i) & ~GnssType::TYPE) + GnssType::RANGE)
              count++;
    }

    return count;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}


/***********************************************/

void GnssSignalBias2SinexBias::run(Config &config)
{
  try
  {
    FileName outNameSinexBias;
    std::vector<FileName> inNameTransmitterInfo;
    std::vector<Data> transmitterBiases, receiverBiases;
    Time timeStart, timeEnd;
    std::string agencyCode, fileAgencyCode, biasMode, method, receiverClockReferenceGnss, description, output, contact, software, hardware, input;
    std::vector<std::string> satelliteClockReferenceObservables, comments;
    UInt sampling = 0;
    UInt intervalLength = 0;

    readConfig(config, "outputfileSinexBias",      outNameSinexBias,      Config::MUSTSET,  "", "");
    readConfig(config, "inputfileTransmitterInfo", inNameTransmitterInfo, Config::MUSTSET,  "", "one file per satellite");
    readConfig(config, "transmitterBiases",        transmitterBiases,     Config::OPTIONAL, "", "one element per satellite");
    readConfig(config, "receiverBiases",           receiverBiases,        Config::OPTIONAL, "", "one element per station");
    readConfig(config, "agencyCode",               agencyCode,            Config::MUSTSET,  "TUG", "identify the agency providing the data");
    readConfig(config, "fileAgencyCode",           fileAgencyCode,        Config::MUSTSET,  "TUG", "identify the agency creating the file");
    readConfig(config, "timeStart",                timeStart,             Config::MUSTSET,  "", "start time of the data");
    readConfig(config, "timeEnd",                  timeEnd,               Config::MUSTSET,  "", "end time of the data ");
    if(readConfigChoice(config, "biasMode", biasMode, Config::MUSTSET, "absolute",  "absolute or relative bias estimates"))
    {
      readConfigChoiceElement(config, "absolute", biasMode, "");
      readConfigChoiceElement(config, "relative", biasMode, "");
      endChoice(config);
    }
    readConfig(config, "observationSampling", sampling,       Config::OPTIONAL, "30", "[seconds]");
    readConfig(config, "intervalLength",      intervalLength, Config::OPTIONAL, "86400", "[seconds] interval for bias parameter representation");
    readConfig(config, "determinationMethod", method,         Config::MUSTSET,   "CO-ESTIMATED_IN_LSA", "determination method used to generate the bias results (see SINEX Bias format description)");
    readConfig(config, "receiverClockReferenceGnss",         receiverClockReferenceGnss,         Config::OPTIONAL, "G", "(G, R, E, C) reference GNSS used for receiver clock estimation");
    readConfig(config, "satelliteClockReferenceObservables", satelliteClockReferenceObservables, Config::OPTIONAL, "G  C1W  C2W", "one per system, reference code observable on first and second frequency (RINEX3 format)");
    readConfig(config, "description",         description,    Config::MUSTSET,   "Graz University of Technology, Austria (TUG/TU Graz)", "organizition gathering/altering the file contents");
    readConfig(config, "contact",             contact,        Config::MUSTSET,   "", "contact name and/or email address");
    readConfig(config, "input",               input,          Config::MUSTSET,   "GNSS observations from the IGS station network", "brief description of the input used to generate this solution");
    readConfig(config, "output",              output,         Config::MUSTSET,   "Estimated signal biases from GNSS processing", "description of the file contents");
    readConfig(config, "software",            software,       Config::MUSTSET,   "GROOPS", "software used to generate the file");
    readConfig(config, "hardware",            hardware,       Config::MUSTSET,   "", "computer hardware on which above software was run");
    readConfig(config, "comment",             comments,       Config::OPTIONAL,  "", "comments in the comment block");
    if(isCreateSchema(config)) return;

    logStatus << "read transmitter info files" << Log::endl;
    for(const auto &fileName : inNameTransmitterInfo)
    {
      GnssStationInfo info;
      readFileGnssStationInfo(fileName, info);
      prn2info[info.markerNumber] = info;
    }

    logStatus << "read bias files" << Log::endl;
    readData(transmitterBiases);
    readData(receiverBiases);

    // create map: GLONASS frequencyNumber -> list of PRNs
    for(const auto &trans : transmitterBiases)
      if(trans.identifier.at(0) == 'R') // GLONASS only
      {
        UInt idRecv = prn2info[trans.identifier].findReceiver(0.5*(timeStart+timeEnd));
        if(idRecv == NULLINDEX)
          throw(Exception("no receiver info found in transmitter info for "+trans.identifier));

        freqNo2prns[std::stoi(prn2info[trans.identifier].receiver.at(idRecv).version)].push_back(trans.identifier);
      }

    logStatus << "write SINEX bias file <" << outNameSinexBias << ">" << Log::endl;
    OutFile file(outNameSinexBias);

    // SINEX header line
    const UInt estimateCount = countBiases(transmitterBiases) + countBiases(receiverBiases);
    agencyCode.resize(3, ' ');
    fileAgencyCode.resize(3, ' ');
    std::transform(biasMode.begin(), biasMode.end(), biasMode.begin(), ::toupper);
    file << "%=BIA " << 1%"%4.2f"s << ' ' << fileAgencyCode << ' ' << Sinex::time2str(System::now(), TRUE) << ' ' << agencyCode << ' '
         << Sinex::time2str(timeStart, TRUE) << ' ' << Sinex::time2str(timeEnd, TRUE) << ' '  << biasMode.at(0) << ' ' << estimateCount%"%08i"s << std::endl;
    file << "*-------------------------------------------------------------------------------" << std::endl;
    file << "* Bias Solution INdependent EXchange Format (Bias-SINEX)" << std::endl;

    // Block: FILE/REFERENCE
    description.resize(std::min(description.size(), UInt(60)), ' ');
    contact.resize(std::min(contact.size(), UInt(60)), ' ');
    input.resize(std::min(input.size(), UInt(60)), ' ');
    output.resize(std::min(output.size(), UInt(60)), ' ');
    software.resize(std::min(software.size(), UInt(60)), ' ');
    hardware.resize(std::min(hardware.size(), UInt(60)), ' ');
    file << "*-------------------------------------------------------------------------------" << std::endl;
    file << "+FILE/REFERENCE" << std::endl;
    file << "*INFO_TYPE_________ INFO________________________________________________________" << std::endl;
    file << " DESCRIPTION        " << description << std::endl;
    file << " CONTACT            " << contact     << std::endl;
    file << " INPUT              " << input       << std::endl;
    file << " OUTPUT             " << output      << std::endl;
    file << " SOFTWARE           " << software    << std::endl;
    file << " HARDWARE           " << hardware    << std::endl;
    file << "-FILE/REFERENCE" << std::endl;

    // Block: FILE/COMMENT
    if(comments.size())
    {
      file << "*-------------------------------------------------------------------------------" << std::endl;
      file << "+FILE/COMMENT" << std::endl;
      for(auto comment : comments)
      {
        comment.resize(std::min(comment.size(), UInt(79)));
        file << ' ' << comment << std::endl;
      }
      file << "-FILE/COMMENT" << std::endl;
    }

    // Block: INPUT/ACKNOWLEDGMENTS

    // Block: BIAS/DESCRIPTION
    method.resize(std::min(method.size(), UInt(39)), ' ');
    file << "*-------------------------------------------------------------------------------" << std::endl;
    file << "+BIAS/DESCRIPTION" << std::endl;
    file << "*KEYWORD________________________________ VALUE(S)_______________________________" << std::endl;
    if(sampling > 0)
      file << " OBSERVATION_SAMPLING                    " << sampling%"% 12i"s << std::endl;
    if(intervalLength > 0)
      file << " PARAMETER_SPACING                       " << intervalLength%"% 12i"s << std::endl;
    file << " DETERMINATION METHOD                    " << method << std::endl;
    file << " BIAS_MODE                               " << biasMode << std::endl;
    file << " TIME_SYSTEM                             G" << std::endl;
    if(!receiverClockReferenceGnss.empty())
    {
      file << " RECEIVER_CLOCK_REFERENCE_GNSS           " << receiverClockReferenceGnss.at(0) << std::endl;
    }
    for(auto sys : satelliteClockReferenceObservables)
    {
      sys.resize(11, ' ');
      file << " SATELLITE_CLOCK_REFERENCE_OBSERVABLES   " << sys << std::endl;
    }
    file << "-BIAS/DESCRIPTION" << std::endl;

    // Block: BIAS/RECEIVER_INFORMATION

    // Block: BIAS/SOLUTION
    file << "*-------------------------------------------------------------------------------" << std::endl;
    file << "+BIAS/SOLUTION" << std::endl;
    file << "*BIAS SVN_ PRN STATION__ OBS1 OBS2 BIAS_START____ BIAS_END______ UNIT __ESTIMATED_VALUE____ _STD_DEV___ __ESTIMATED_SLOPE____ _STD_DEV___" << std::endl;
    for(const auto &data : transmitterBiases)
    {
      writeData(file, timeStart, timeEnd, data, /*isStation*/FALSE);
    }
    for(const auto &data : receiverBiases)
      writeData(file, timeStart, timeEnd, data, /*isStation*/TRUE);
    file << "-BIAS/SOLUTION" << std::endl;

    file << "%=ENDBIA" << std::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
