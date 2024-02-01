/***********************************************/
/**
* @file sinexMetadata2GlonassFrequencyNumber.cpp
*
* @brief Create glonassPrnSvn2FrequencyNumber matrix from IGS SINEX metadata file.
*
* @author Torsten Mayer-Guerr
* @date 2023-12-20
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Create \configFile{outputfileMatrixPrn2FrequencyNumber}{matrix} matrix
from \href{https://www.igs.org/mgex/metadata/#metadata}{IGS SINEX metadata format}
with the columns: GLONASS PRN, SVN, mjdStart, mjdEnd, frequencyNumber.

See also \program{RinexObservation2GnssReceiver}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "inputOutput/file.h"
#include "inputOutput/fileSinex.h"

/***** CLASS ***********************************/

/** @brief Create glonassPrnSvn2FrequencyNumber matrix from IGS SINEX metadata file.
* @ingroup programsConversionGroup */
class SinexMetadata2GlonassFrequencyNumber
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(SinexMetadata2GlonassFrequencyNumber, SINGLEPROCESS, "Create glonassPrnSvn2FrequencyNumber matrix from IGS SINEX metadata file.", Conversion, Gnss)

/***********************************************/

void SinexMetadata2GlonassFrequencyNumber::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNamePrn2FrequencyNumber, fileNameSinexMetadata;

    readConfig(config, "outputfileMatrixPrn2FrequencyNumber", fileNamePrn2FrequencyNumber,  Config::MUSTSET, "{groopsDataDir}/gnss/transmitter/glonassPrnSvn2FrequencyNumber.txt", "GROOPS matrix with columns: GLONASS PRN, SVN, mjdStart, mjdEnd, frequencyNumber");
    readConfig(config, "inputfileSinexMetadata",              fileNameSinexMetadata,        Config::MUSTSET, "igs_satellite_metadata.snx", "IGS SINEX metadata file");
    if(isCreateSchema(config)) return;

    Sinex sinex;
    readFileSinex(fileNameSinexMetadata, sinex);
    std::vector<std::string> freqLines = sinex.findBlock("SATELLITE/FREQUENCY_CHANNEL")->lines;
    std::vector<std::string> prnLines  = sinex.findBlock("SATELLITE/PRN")->lines;

    std::vector<Vector> lines;
    for(const auto &line : freqLines)
    {
      // *         1         2         3         4         5         6         7         8
      // *12345678901234567890123456789012345678901234567890123456789012345678901234567890
      // *SVN_ Valid_From____ Valid_To______ chn Comment________________________________
      //  R701 2003:344:00000 2009:258:86399   1 [FC10]
      std::string  svn = line.substr(1, 4);
      Time   timeStart = Sinex::str2time(line,  6, FALSE/*zeroIsMaxTime*/, TRUE/*fourDigitYear*/);
      Time   timeEnd   = Sinex::str2time(line, 21, TRUE /*zeroIsMaxTime*/, TRUE/*fourDigitYear*/);
      Double chn       = String::toDouble(line.substr(36, 3));

      for(const auto &line : prnLines)
        if(line.substr(1, 4) == svn)
        {
          // *         1         2         3         4         5         6         7         8
          // *12345678901234567890123456789012345678901234567890123456789012345678901234567890
          // *SVN_ Valid_From____ Valid_To______ PRN Comment_________________________________
          //  G001 1978:053:00000 1985:199:00000 G04
          std::string prn = line.substr(36, 3);
          Time tStart  = Sinex::str2time(line,  6, FALSE/*zeroIsMaxTime*/, TRUE/*fourDigitYear*/);
          Time tEnd    = Sinex::str2time(line, 21, TRUE /*zeroIsMaxTime*/, TRUE/*fourDigitYear*/);
          if((tStart < timeEnd) && (tEnd > timeStart))
            lines.push_back(Vector({String::toDouble(prn.substr(1)), String::toDouble(svn.substr(1)), std::max(timeStart, tStart).mjd(), std::min(timeEnd, tEnd).mjd(), chn}));
        }
    }

    std::sort(lines.begin(), lines.end(), [](auto &x, auto &y) {if(x(0) != y(0)) return x(0)<y(0); return x(2)<y(2);});

    logStatus<<"write frequency number file <"<<fileNamePrn2FrequencyNumber<<">"<<Log::endl;
    OutFile file(fileNamePrn2FrequencyNumber);
    file<<"# PRN SVN  MJD start          MJD end           chn"<<std::endl;
    file<<"# -------------------------------------------------"<<std::endl;
    for(const auto &line : lines)
      file<<line(0)%" %02i"s<<line(1)%"  %03i"s<<line(2)%" %18.10f"s<<line(3)%" %18.10f"s<<line(4)%" %4i"s<<std::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
