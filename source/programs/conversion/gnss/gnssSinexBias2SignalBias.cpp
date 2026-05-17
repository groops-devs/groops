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
Time-variable biases (slopes) are ignored. Combined \verb|X| signals are transformed to actually transmitted signals.

If entries exist for different time periods, you must select one using \config{time}.

See also \program{GnssSignalBias2SinexBias}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "inputOutput/fileSinex.h"
#include "files/fileMatrix.h"
#include "files/fileGnssSignalBias.h"
#include "gnss/gnssReceiver.h"

/***** CLASS ***********************************/

/** @brief Converts GNSS signal biases from SINEX Bias format to GnssSignalBias format.
* @ingroup programsGroup */
class GnssSinexBias2SignalBias
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GnssSinexBias2SignalBias, SINGLEPROCESS, "Converts GNSS signal biases from SINEX Bias format to GnssSignalBias format.", Conversion, Gnss)

/***********************************************/

void GnssSinexBias2SignalBias::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameSignalBias, fileNameSinexBias, fileNamePrnSvn2FreqNo;
    std::vector<std::string> identifiers;
    Time time;

    readConfig(config, "outputfileSignalBias",         fileNameSignalBias,    Config::MUSTSET,  "", "prn is appended to file name");
    readConfig(config, "inputfileSinexBias",           fileNameSinexBias,     Config::MUSTSET,  "", "");
    readConfig(config, "inputfilePrn2FrequencyNumber", fileNamePrnSvn2FreqNo, Config::OPTIONAL, "{groopsDataDir}/gnss/transmitter/glonassPrnSvn2FrequencyNumber.txt", "matrix with columns: GLONASS PRN, SVN, mjdStart, mjdEnd, frequencyNumber");
    readConfig(config, "identifier",                   identifiers,           Config::OPTIONAL, "", "(empty = all) satellite PRN, e.g. G23 or E05");
    readConfig(config, "time",                         time,                  Config::OPTIONAL, "", "if entries exist for different time periods");
    if(isCreateSchema(config)) return;

    // ==================================================

    Matrix prnSvn2FreqNo;
    if(!fileNamePrnSvn2FreqNo.empty())
    {
      logStatus<<"read GLONASS PRN/SVN to frequency number matrix <"<<fileNamePrnSvn2FreqNo<<">"<< Log::endl;
      readFileMatrix(fileNamePrnSvn2FreqNo, prnSvn2FreqNo);
    }

    // ==================================================

    logStatus<<"read SINEX Bias file <"<<fileNameSinexBias<<">"<<Log::endl;
    logInfo<<"  (only satellite specific OSB (observable-specific signal bias) are considered)"<<Log::endl;
    Sinex sinex;
    readFileSinex(fileNameSinexBias, sinex);

    std::map<std::string, GnssSignalBias> prn2signalBias;
    for(const auto &line : sinex.findBlock("BIAS/SOLUTION")->lines)
    {
      // *         1         2         3         4         5         6         7         8        9         10        11        12        13
      // *123456789012345678901234567890123456789012345678901234567890123456789012345678902345678901234567890123456789012345678901234567890123456
      // *BIAS SVN_ PRN STATION__ OBS1 OBS2 BIAS_START____ BIAS_END______ UNIT __ESTIMATED_VALUE____ _STD_DEV___ __ESTIMATED_SLOPE____ _STD_DEV__
      std::string biasType  = String::trim(line.substr(1,3));
      std::string svn       = String::trim(line.substr(6,4));
      std::string prn       = String::trim(line.substr(11,3));
      std::string station   = String::trim(line.substr(15,9));
      GnssType    gnssType  = GnssType(line.substr(25,3)+prn);
      Time        timeStart = Sinex::str2time(line, 35, FALSE, TRUE);
      Time        timeEnd   = Sinex::str2time(line, 50, TRUE,  TRUE);
      Time        timeMid   = 0.5*(timeStart+timeEnd);
      std::string unit      = String::trim(line.substr(65,4));
      Double      value     = String::toDouble(line.substr(70,21));
   // Double      slope     = (line.size() >= 126) ? String::toDouble(line.substr(105,21)) : 0.;

      // satellite specific OSB only
      if((biasType != "OSB") || prn.empty() || !station.empty())
        continue;
      if(identifiers.size() && std::find(identifiers.begin(), identifiers.end(), prn) == identifiers.end())
        continue;
      if((time != Time()) && !((timeStart <= time) && (time <= timeEnd)))
        continue;

      // add GLONASS frequency number (find PRN, SVN, time interval)
      if(gnssType == GnssType::GLONASS)
        for(UInt i=0; i<prnSvn2FreqNo.rows(); i++)
          if((prn == prnSvn2FreqNo(i,0)%"R%02i"s) && (svn == prnSvn2FreqNo(i,1)%"R%03i"s) &&
             (mjd2time(prnSvn2FreqNo(i,2)) <= timeMid) && (timeMid <= mjd2time(prnSvn2FreqNo(i,3))))
            gnssType.setFrequencyNumber(static_cast<Int>(prnSvn2FreqNo(i,4)));

      if(unit == "ns")
        value *= 1e-9*LIGHT_VELOCITY;
      else if(unit == "cyc")
        value *= gnssType.wavelength();
      else
        throw(Exception("unknown unit: "+unit));

      if(gnssType == GnssType::PHASE)
        gnssType = gnssType & ~GnssType::ATTRIBUTE;

      UInt idx;
      if(!gnssType.isInList(prn2signalBias[prn].types, idx))
      {
        prn2signalBias[prn].types.push_back(gnssType);
        prn2signalBias[prn].biases.push_back(value);
      }
      else if(value != prn2signalBias[prn].biases.at(idx))
        logWarning<<prn<<": different biases for the same type "<<gnssType.str()<<" -> only the first one is used"<<Log::endl;
    }

    if(prn2signalBias.empty())
    {
      logWarning<<"no satellites found"<<Log::endl;
      return;
    }

    // ==================================================

    logStatus<<"write signal bias files <"<<fileNameSignalBias<<">"<<Log::endl;
    for(const auto &prn : prn2signalBias)
    {
      GnssSignalBias signalBias;

      // Transform combined X signals to actually transmitted signals.
      Matrix T;
      GnssReceiver::signalCompositionDefault(prn.second.types, signalBias.types, T);
      T = pseudoInverse(T);
      signalBias.biases = Vector(T * Vector(prn.second.biases));

      writeFileGnssSignalBias(fileNameSignalBias.appendBaseName("."+prn.first), signalBias);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
