/***********************************************/
/**
* @file swarm2Starcamera.cpp
*
* @brief read SWARM attitude data.
*
* @author Norbert Zehentner
* @date 2014-05-23
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program reads SWARM star camera data given in the cdf format
and before converted to an ascii file using the program \verb|cdfexport|
provided by the Goddard Space Flight Center (\url{http://cdf.gsfc.nasa.gov/}).
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "files/fileInstrument.h"
#include "classes/earthRotation/earthRotation.h"

/***** CLASS ***********************************/

/** @brief read SWARM star camera data.
* @ingroup programsConversionGroup */
class Swarm2Starcamera
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Swarm2Starcamera, SINGLEPROCESS, "read SWARM star camera data", Conversion, Instrument)

/***********************************************/

void Swarm2Starcamera::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName              fileNameSca;
    std::vector<FileName> fileNamesIn;
    EarthRotationPtr      earthRotation;

    readConfig(config, "outputfileStarCamera", fileNameSca,   Config::MUSTSET,  "", "");
    readConfig(config, "earthRotation",        earthRotation, Config::MUSTSET,  "", "");
    readConfig(config, "inputfile",            fileNamesIn,   Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    logStatus<<"read input files"<<Log::endl;
    StarCameraArc arcSca;
    for(const auto &fileName : fileNamesIn)
    {
      try
      {
        logStatus<<"read file <"<<fileName<<">"<<Log::endl;
        InFile file(fileName);
        file.exceptions(std::ios::badbit|std::ios::failbit);

        // find begin of data
        for(;;)
        {
          std::string line;
          std::getline(file, line);
          if(line.find("Timestamp") != std::string::npos)
            break;
        }

        // read data
        for(;;)
        {
          std::string line;
          try
          {
            std::getline(file, line);
            if(line.empty())
              continue;
            }
            catch(std::exception &/*e*/)
            {
              break;
            }

          const std::vector<std::string> monthStr = {"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov","Dec"};
          UInt   day     = String::toInt(line.substr(0, 2));
          UInt   month   = std::distance(monthStr.begin(), std::find(monthStr.begin(), monthStr.end(), line.substr(3,3))) + 1;
          UInt   year    = String::toInt(line.substr(7, 4));
          UInt   hour    = String::toInt(line.substr(12, 2));
          UInt   minute  = String::toInt(line.substr(15, 2));
          Double seconds = String::toDouble(line.substr(18, 6));

          Vector q(4);
          q(1) = String::toDouble(line.substr(25, 12));
          std::getline(file, line);
          q(2) = String::toDouble(line.substr(25, 12));
          std::getline(file, line);
          q(3) = String::toDouble(line.substr(25, 12));
          std::getline(file, line);
          q(0) = String::toDouble(line.substr(25, 12));

          if(std::fabs(norm(q)-1) > 1e-5)
          {
            logWarning<<"strange norm = "<<norm(q)<<Log::endl;
            continue;
          }

          StarCameraEpoch epoch;
          epoch.time   = timeUTC2GPS(date2time(year, month, day, hour, minute, seconds));
          epoch.rotary = inverse(Rotary3d(q) * earthRotation->rotaryMatrix(epoch.time));
          arcSca.push_back(epoch);
        }
      }
      catch(std::exception &e)
      {
        logError<<e.what()<<": continue..."<<Log::endl;
      }
    } // for(idFile)

    // write data
    // ----------
    logStatus<<"write star camera data to file <"<<fileNameSca<<">"<<Log::endl;
    InstrumentFile::write(fileNameSca, arcSca);
    Arc::printStatistics(arcSca);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
