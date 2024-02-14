/***********************************************/
/**
* @file goceXml2Orbit.cpp
*
* @brief Read ESA XML GOCE Data.
*
* @author Torsten Mayer-Guerr
* @date 2010-10-19
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Read ESA XML GOCE Data.
)";

/***********************************************/

#include "programs/program.h"
#include "parser/xml.h"
#include "files/fileInstrument.h"
#include "classes/earthRotation/earthRotation.h"

/***** CLASS ***********************************/

/** @brief Read ESA XML GOCE Data.
* @ingroup programsConversionGroup */
class GoceXml2Orbit
{
  void readFileGoceOrbit(const FileName &fileName, OrbitArc &arc);

public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GoceXml2Orbit, SINGLEPROCESS, "read ESA XML GOCE Data", Conversion, Orbit, Instrument)

/***********************************************/

void GoceXml2Orbit::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName outName;
    std::vector<FileName> fileName;
    EarthRotationPtr earthRotation;

    readConfig(config, "outputfileOrbit", outName,       Config::MUSTSET,  "", "");
    readConfig(config, "earthRotation",   earthRotation, Config::OPTIONAL, "file", "rotation from TRF to CRF");
    readConfig(config, "inputfile",       fileName,      Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    logStatus<<"read input files"<<Log::endl;
    OrbitArc orbit;
    for(UInt i=0; i<fileName.size(); i++)
    {
      logStatus<<"read file <"<<fileName.at(i)<<">"<<Log::endl;
      readFileGoceOrbit(fileName.at(i), orbit);
    }

    if(earthRotation)
    {
      logStatus<<"rotation from TRF to CRF"<<Log::endl;
      Single::forEach(orbit.size(), [&](UInt i)
      {
        const Rotary3d rotation = inverse(earthRotation->rotaryMatrix(orbit.at(i).time));
        const Vector3d omega    = earthRotation->rotaryAxis(orbit.at(i).time);
        orbit.at(i).position = rotation.rotate(orbit.at(i).position);
        orbit.at(i).velocity = rotation.rotate(orbit.at(i).velocity) + crossProduct(omega, orbit.at(i).position);
      });
    }

    logStatus<<"write orbit data to file <"<<outName<<">"<<Log::endl;
    InstrumentFile::write(outName, orbit);
    Arc::printStatistics(orbit);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GoceXml2Orbit::readFileGoceOrbit(const FileName &fileName, OrbitArc &arc)
{
  try
  {
    OrbitEpoch epoch;

    InFile     file(fileName);
    XmlNodePtr rootNode = XmlNode::read(file);
    XmlNodePtr sstNode  = getChild(rootNode, "SST_PSO_2",  TRUE);
    sstNode  = getChild(sstNode, "SST_PRD_2",  TRUE);
    sstNode  = getChild(sstNode, "List_of_SP3c_Records", TRUE);

    UInt epochCount = childCount(sstNode, "SP3c_Record", TRUE);
    for(UInt i=0; i<epochCount; i++)
    {
      XmlNodePtr epochNode = getChild(sstNode, "SP3c_Record", TRUE);

      UInt   year, month, day, hour, min;
      Double sec = 0;
      XmlNodePtr timeNode  = getChild(epochNode, "Time_Information", TRUE);
      timeNode  = getChild(timeNode,  "GPS_Time",  TRUE);
      timeNode  = getChild(timeNode,  "Start",     TRUE);
      timeNode  = getChild(timeNode,  "Gregorian", TRUE);
      readXml(timeNode, "Year",         year,  TRUE);
      readXml(timeNode, "Month",        month, TRUE);
      readXml(timeNode, "Day_of_Month", day,   TRUE);
      readXml(timeNode, "Hour",         hour,  TRUE);
      readXml(timeNode, "Minute",       min,   TRUE);
      readXml(timeNode, "Second",       sec,   TRUE);
      epoch.time = date2time(year, month, day, hour, min, sec);

      XmlNodePtr satNode = getChild(getChild(epochNode, "List_of_Satellite_IDs", TRUE), "L15", TRUE);
      XmlNodePtr posNode = getChild(satNode, "Position", TRUE);
      readXml(posNode, "X", epoch.position.x(),  TRUE);
      readXml(posNode, "Y", epoch.position.y(),  TRUE);
      readXml(posNode, "Z", epoch.position.z(),  TRUE);
      epoch.position *= 1e3;
      XmlNodePtr velNode = getChild(satNode, "Velocity", TRUE);
      readXml(velNode, "X", epoch.velocity.x(),  TRUE);
      readXml(velNode, "Y", epoch.velocity.y(),  TRUE);
      readXml(velNode, "Z", epoch.velocity.z(),  TRUE);
      epoch.velocity *= 0.1;
      arc.push_back(epoch);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

// void GoceXml2Orbit::readFileGoceOrbit(const FileName &fileName, OrbitArc &arc)
// {
//   try
//   {
//     OrbitEpoch epoch;
//
//     XmlNodePtr rootNode = XmlNode::readFile(fileName);
//     XmlNodePtr dataNode = getChild(rootNode, "Data_Block", TRUE);
//     XmlNodePtr sstNode  = getChild(dataNode, "SST_PVT_DS", TRUE);
//
//     UInt epochCount = childCount(sstNode, "SST_PVT_1i", TRUE);
//     for(UInt i=0; i<epochCount; i++)
//     {
//       XmlNodePtr epochNode = getChild(sstNode, "SST_PVT_1i", TRUE);
//       Double t = 0;
//       readXml(epochNode, "Tt_GPS", t, TRUE);
//       epoch.time = seconds2time(t) + date2time(1980,1,6);
//
//       XmlNodePtr posNode = getChild(epochNode, "SST_Pos", TRUE);
//       std::string posStr;
//       readXml(posNode, "Position", posStr, TRUE);
//       std::stringstream ss(posStr);
//       ss>>epoch.position.x()>>epoch.position.y()>>epoch.position.z();
//
//       std::string velStr;
//       readXml(epochNode, "SST_Vel", velStr, TRUE);
//       std::stringstream ss2(velStr);
//       ss2>>epoch.velocity.x()>>epoch.velocity.y()>>epoch.velocity.z();
//
//       arc.push_back(epoch);
//     }
//   }
//   catch(std::exception &e)
//   {
//     GROOPS_RETHROW(e)
//   }
// }

/***********************************************/
