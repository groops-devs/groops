/***********************************************/
/**
* @file jason2Starcamera.cpp
*
* @brief read Jason star camera data.
*
* @author Norbert Zehentner
* @date 2014-10-16
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program reads in Jason star camera data given in a special  format.
Files available at: \url{cddis.gsfc.nasa.gov/pub/doris/ancillary/quaternions/ja2/}.
A description of the format can be found under:
\url{ftp://ftp.ids-doris.org/pub/ids/ancillary/quaternions/jason1_2_quaternion_solar_panel.pdf}
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "inputOutput/file.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief read Jason star camera data.
* @ingroup programsConversionGroup */
class Jason2Starcamera
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Jason2Starcamera, SINGLEPROCESS, "read Jason star camera data", Conversion, Instrument)

/***********************************************/

void Jason2Starcamera::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName              outNameStar;
    std::vector<FileName> fileNamesIn;
    UInt                  jasonNumber;

    readConfig(config, "outputfileStarCamera", outNameStar, Config::MUSTSET,  "", "");
    readConfig(config, "jasonNumber",          jasonNumber, Config::MUSTSET, "2", "Jason number (different file format), 1 for Sentinel");
    readConfig(config, "inputfile",            fileNamesIn, Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    logStatus<<"read input files"<<Log::endl;
    StarCameraArc arc;
    for(const auto &fileName : fileNamesIn)
    {
      try
      {
        logStatus<<"read file <"<<fileName<<">"<<Log::endl;
        InFile file(fileName);

        std::string line;
        while(std::getline(file, line))
        {
          if(line.empty())
            continue;

          std::string lineID1 = line.substr(0,1);
          if(lineID1 == "#" || lineID1 == "Q" || (line.length()<25))
            continue;

          const UInt   year   = String::toInt(line.substr(0, 4));
          const UInt   month  = String::toInt(line.substr(5, 2));
          const UInt   day    = String::toInt(line.substr(8, 2));
          const UInt   hour   = String::toInt(line.substr(11, 2));
          const UInt   minute = String::toInt(line.substr(14, 2));
          const Double second = String::toDouble(line.substr(17, 6));
          if((year < 1990) || (year > 2030))
            continue;

          std::stringstream ss(line.substr(23));
          ss.exceptions(std::ios::badbit | std::ios::failbit);

          Vector q(4);
          Double tmp;
          if(jasonNumber == 1)
            ss>>q(0)>>q(1)>>q(2)>>q(3);
          else
            ss>>tmp>>q(0)>>tmp>>tmp>>q(1)>>tmp>>tmp>>q(2)>>tmp>>tmp>>q(3)>>tmp;

          if(fabs(norm(q)-1) > 1e-5)
          {
            logWarning<<timeUTC2GPS(date2time(year, month, day, hour, minute, second)).dateTimeStr()<<" strange quaternion norm = "<<norm(q)<<Log::endl;
            continue;
          }

          StarCameraEpoch epoch;
          epoch.time   = timeUTC2GPS(date2time(year, month, day, hour, minute, second));
          epoch.rotary = Rotary3d(q);
          arc.push_back(epoch);
        }
      }
      catch(std::exception &e)
      {
        logError<<e.what()<<": continue..."<<Log::endl;
      }
    }

    logStatus<<"write star camera data to file <"<<outNameStar<<">"<<Log::endl;
    InstrumentFile::write(outNameStar, arc);
    Arc::printStatistics(arc);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

