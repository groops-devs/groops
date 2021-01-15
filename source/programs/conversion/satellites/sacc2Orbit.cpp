/***********************************************/
/**
* @file sacc2Orbit.cpp
*
* @brief read SACC orbit data.
*
* @author Norbert Zehentner
* @date 2014-08-22
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program reads in SACC orbit data.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief read SACC orbit data.
* @ingroup programsConversionGroup */
class Sacc2Orbit
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Sacc2Orbit, SINGLEPROCESS, "read SACC orbit data", Conversion, Orbit, Instrument)

/***********************************************/

void Sacc2Orbit::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName              fileNameOrbit;
    std::vector<FileName> fileNamesIn;

    readConfig(config, "outputfileOrbit", fileNameOrbit, Config::MUSTSET,  "", "");
    readConfig(config, "inputfile",       fileNamesIn,   Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    logStatus<<"read input files"<<Log::endl;
    OrbitArc arc;
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
          std::stringstream ss(line);
          ss.exceptions(std::ios::badbit | std::ios::failbit);

          UInt       id, year, month, day, hour, minute;
          Double     second;
          OrbitEpoch epoch;

          ss>>id>>year>>month>>day>>hour>>minute>>second;
          ss>>epoch.position.x()>>epoch.position.y()>>epoch.position.z();
          ss>>epoch.velocity.x()>>epoch.velocity.y()>>epoch.velocity.z();

          epoch.time = timeUTC2GPS(date2time(year, month, day, hour, minute, second));
          epoch.position *= 1000; // km => m
          epoch.velocity *= 1000; // km/s => m/s
          arc.push_back(epoch);
        } // while(getline)
      }
      catch(std::exception &e)
      {
        logError<<e.what()<<": continue..."<<Log::endl;
      }
    } // for(idFile)

    if(!fileNameOrbit.empty())
    {
      logInfo<<"write data to <"<<fileNameOrbit<<">"<<Log::endl;
      InstrumentFile::write(fileNameOrbit, arc);
      Arc::printStatistics(arc);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

