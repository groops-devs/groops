/***********************************************/
/**
* @file sentinel2orbit.cpp
*
* @brief Read Sentinel orbits from XML format.
*
* @author Barbara Suesser-Rechberger
* @date 2019-11-25
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Read Sentinel orbits from XML format.
)";

/***********************************************/

#include "programs/program.h"
#include "parser/xml.h"
#include "files/fileInstrument.h"
#include "classes/earthRotation/earthRotation.h"


/***** CLASS ***********************************/

/** @brief Read XML Sentinel Data.
* @ingroup programsConversionGroup */
class SentinelXml2Orbit
{
public:
  void run(Config& config);
};

GROOPS_REGISTER_PROGRAM(SentinelXml2Orbit, SINGLEPROCESS, "read XML Sentinel Data", Conversion, Orbit, Instrument)

/***********************************************/

void SentinelXml2Orbit::run(Config& config)
{
  try
  {
    FileName              fileNameOrbit;
    std::vector<FileName> fileNamesIn;
    EarthRotationPtr      earthRotation;

    readConfig(config, "outputfileOrbit", fileNameOrbit, Config::MUSTSET,  "", "");
    readConfig(config, "earthRotation",   earthRotation, Config::MUSTSET,  "", "");
    readConfig(config, "inputfile",       fileNamesIn,   Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    logStatus<<"read input files"<<Log::endl;
    OrbitArc orbit;
    for(const auto &fileName : fileNamesIn)
    {
      try
      {
        logStatus<<"read file <"<<fileName<<">"<<Log::endl;
        InFile file(fileName);
        XmlNodePtr rootNode    = XmlNode::read(file);
        XmlNodePtr dataNode    = getChild(rootNode, "Data_Block",   TRUE);
        XmlNodePtr osvListNode = getChild(dataNode, "List_of_OSVs", TRUE);

        const UInt epochCount = childCount(osvListNode, "OSV", TRUE);
        for(UInt i=0; i<epochCount; i++)
        {
          XmlNodePtr epochNode = getChild(osvListNode, "OSV", TRUE);

          // Read UTC string and convert it into modified julian date and then into GPS time
          std::string timeUTC;
          Char   tmp;
          UInt   year, month, day, hour, minute;
          Double second;
          readXml(epochNode, "UTC", timeUTC, TRUE);
          std::stringstream ss(timeUTC.substr(4));
          ss>>year>>tmp>>month>>tmp>>day>>tmp>>hour>>tmp>>minute>>tmp>>second;

          OrbitEpoch epoch;
          epoch.time = timeUTC2GPS(date2time(year, month, day, hour, minute, second));

          // position coordinates
          readXml(epochNode, "X", epoch.position.x(), TRUE);
          readXml(epochNode, "Y", epoch.position.y(), TRUE);
          readXml(epochNode, "Z", epoch.position.z(), TRUE);

          // Velocity coordinates
          readXml(epochNode, "VX", epoch.velocity.x(), TRUE);
          readXml(epochNode, "VY", epoch.velocity.y(), TRUE);
          readXml(epochNode, "VZ", epoch.velocity.z(), TRUE);

          orbit.push_back(epoch);
        }
      }
      catch(std::exception &e)
      {
        logError<<e.what()<<": continue..."<<Log::endl;
      }
    }

    // Rotation TRF -> CRF
    // -------------------
    if(earthRotation)
    {
      logStatus<<"rotation from TRF to CRF"<<Log::endl;
      logTimerStart;
      for(UInt i=0; i<orbit.size(); i++)
      {
        logTimerLoop(i, orbit.size());
        const Rotary3d rotation = inverse(earthRotation->rotaryMatrix(orbit.at(i).time));
        orbit.at(i).position    = rotation.rotate(orbit.at(i).position);
        if(orbit.at(i).velocity.r() > 0)
          orbit.at(i).velocity = rotation.rotate(orbit.at(i).velocity) + crossProduct(earthRotation->rotaryAxis(orbit.at(i).time), orbit.at(i).position);
      }
      logTimerLoopEnd(orbit.size());
    }

    if(!fileNameOrbit.empty())
    {
      logInfo<<"write data to <"<<fileNameOrbit<<">"<<Log::endl;
      InstrumentFile::write(fileNameOrbit, orbit);
      Arc::printStatistics(orbit);
    }
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
