/***********************************************/
/**
* @file champ2Orbit.cpp
*
* @brief read champ PSO orbits from the special CHORB format
*
* @author Norbert Zehentner
* @date 2011-10-04
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program reads in CHAMP precise science orbits in the special CHORB format.
A description of the format can be found under: \url{http://op.gfz-potsdam.de/champ/docs_CHAMP/CH-GFZ-FD-002.pdf}
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "files/fileInstrument.h"
#include "classes/earthRotation/earthRotation.h"

/***** CLASS ***********************************/

/** @brief read champ PSO orbits from the special CHORB format
* @ingroup programsConversionGroup */
class Champ2Orbit
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Champ2Orbit, SINGLEPROCESS, "read champ PSO orbits from the special CHORB format", Conversion, Orbit, Instrument)

/***********************************************/

void Champ2Orbit::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName              fileNameOrbit;
    std::vector<FileName> fileNamesInput;
    EarthRotationPtr      earthRotation;

    readConfig(config, "outputfileOrbit", fileNameOrbit,  Config::MUSTSET,  "",     "");
    readConfig(config, "earthRotation",   earthRotation,  Config::OPTIONAL, "file", "");
    readConfig(config, "inputfile",       fileNamesInput, Config::MUSTSET,  "",     "");
    if(isCreateSchema(config)) return;

    OrbitArc arc;
    for(auto &fileName : fileNamesInput)
    {
      logStatus<<"read file <"<<fileName<<">"<<Log::endl;
      InFile file(fileName);
      std::string line;
      while(std::getline(file, line))
        if(String::startsWith(line, "ORBIT"))
          break;

      while(std::getline(file, line))
      {
        OrbitEpoch epoch;
        UInt   day = String::toInt(line.substr(0, 6));
        Double ms  = String::toDouble(line.substr(6, 11));
        epoch.time         = timeTT2GPS(mjd2time(J2000+0.1*day) + seconds2time(1e-6*ms));
        epoch.position.x() = 1e-3*String::toDouble(line.substr(17, 12));
        epoch.position.y() = 1e-3*String::toDouble(line.substr(29, 12));
        epoch.position.z() = 1e-3*String::toDouble(line.substr(41, 12));
        epoch.velocity.x() = 1e-7*String::toDouble(line.substr(53, 12));
        epoch.velocity.y() = 1e-7*String::toDouble(line.substr(65, 12));
        epoch.velocity.z() = 1e-7*String::toDouble(line.substr(77, 12));
        arc.push_back(epoch);
      }
    }

    if(earthRotation)
    {
      logStatus<<"rotation from TRF to CRF"<<Log::endl;
      Single::forEach(arc.size(), [&](UInt i)
      {
        const Rotary3d rotation = inverse(earthRotation->rotaryMatrix(arc.at(i).time));
        arc.at(i).position = rotation.rotate(arc.at(i).position);
        if(arc.at(i).velocity.r() > 0)
          arc.at(i).velocity = rotation.rotate(arc.at(i).velocity) + crossProduct(earthRotation->rotaryAxis(arc.at(i).time), arc.at(i).position);
      });
    }

    if(!fileNameOrbit.empty())
    {
      logStatus<<"write orbit data to file <"<<fileNameOrbit<<">"<<Log::endl;
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
