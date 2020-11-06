/***********************************************/
/**
* @file settings.cpp
*
* @brief Read/write constants.
*
* @author Torsten Mayer-Guerr
* @date 2010-03-01
*
*/
/***********************************************/

#include "base/importStd.h"
#include "base/constants.h"
#include "base/time.h"
#include "parser/xml.h"
#include "inputOutput/fileName.h"
#include "inputOutput/file.h"
#include "settings.h"

/***** FUNCTIONS *******************************/

void readFileSettings(const FileName &fileName)
{
  try
  {
    InFile stream(fileName);
    XmlNodePtr xmlNode = XmlNode::read(stream);

    readXml(xmlNode, "LIGHT_VELOCITY", LIGHT_VELOCITY);
    readXml(xmlNode, "DEFAULT_GM",     STRING_DEFAULT_GM);
    readXml(xmlNode, "DEFAULT_R",      STRING_DEFAULT_R);
    readXml(xmlNode, "GRS80_a",        STRING_DEFAULT_GRS80_a);
    readXml(xmlNode, "GRS80_f",        STRING_DEFAULT_GRS80_f);
    DEFAULT_GM      = std::atof(STRING_DEFAULT_GM.c_str());
    DEFAULT_R       = std::atof(STRING_DEFAULT_R.c_str());
    DEFAULT_GRS80_a = std::atof(STRING_DEFAULT_GRS80_a.c_str());
    DEFAULT_GRS80_f = std::atof(STRING_DEFAULT_GRS80_f.c_str());
    readXml(xmlNode, "GRAVITATIONALCONSTANT", GRAVITATIONALCONSTANT);
    readXml(xmlNode, "R_Earth",       R_Earth);
    readXml(xmlNode, "R_Moon",        R_Moon);
    readXml(xmlNode, "GM_Earth",      GM_Earth);
    readXml(xmlNode, "GM_Sun",        GM_Sun);
    readXml(xmlNode, "GM_Moon",       GM_Moon);
    readXml(xmlNode, "GM_MERCURY",    GM_MERCURY);
    readXml(xmlNode, "GM_VENUS",      GM_VENUS);
    readXml(xmlNode, "GM_MARS",       GM_MARS);
    readXml(xmlNode, "GM_JUPITER",    GM_JUPITER);
    readXml(xmlNode, "GM_SATURN",     GM_SATURN);
    readXml(xmlNode, "TIME_EPSILON",  TIME_EPSILON);
    readXml(xmlNode, "DELTA_TAI_GPS", DELTA_TAI_GPS);
    readXml(xmlNode, "DELTA_TT_GPS",  DELTA_TT_GPS);
    readXml(xmlNode, "J2000",         STRING_J2000);
    J2000 = atof(STRING_J2000.c_str());

    // read leap seconds
    UInt count = xmlNode->getChildCount("leapSecond");
    MJD_UTC_GPS.resize(count);
    DELTA_UTC_GPS.resize(count);
    for(UInt i=0; i<count; i++)
    {
      XmlNodePtr xmlNode2 = getChild(xmlNode, "leapSecond");
      readXml(xmlNode2, "MJD",            MJD_UTC_GPS.at(i));
      readXml(xmlNode2, "DELTA_UTC_GPS",  DELTA_UTC_GPS.at(i));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void writeFileSettings(const FileName &fileName)
{
  try
  {
    XmlNodePtr xmlNode = XmlNode::create("groops");

    writeXml(xmlNode, "LIGHT_VELOCITY", LIGHT_VELOCITY);
    writeXml(xmlNode, "DEFAULT_GM",     STRING_DEFAULT_GM);
    writeXml(xmlNode, "DEFAULT_R",      STRING_DEFAULT_R);
    writeXml(xmlNode, "GRS80_a",        STRING_DEFAULT_GRS80_a);
    writeXml(xmlNode, "GRS80_f",        STRING_DEFAULT_GRS80_f);
    writeXml(xmlNode, "GRAVITATIONALCONSTANT", GRAVITATIONALCONSTANT);
    writeXml(xmlNode, "R_Earth",        R_Earth);
    writeXml(xmlNode, "R_Moon",         R_Moon);
    writeXml(xmlNode, "GM_Earth",       GM_Earth);
    writeXml(xmlNode, "GM_Sun",         GM_Sun);
    writeXml(xmlNode, "GM_Moon",        GM_Moon);
    writeXml(xmlNode, "GM_MERCURY",     GM_MERCURY);
    writeXml(xmlNode, "GM_VENUS",       GM_VENUS);
    writeXml(xmlNode, "GM_MARS",        GM_MARS);
    writeXml(xmlNode, "GM_JUPITER",     GM_JUPITER);
    writeXml(xmlNode, "GM_SATURN",      GM_SATURN);
    writeXml(xmlNode, "TIME_EPSILON",   TIME_EPSILON);
    writeXml(xmlNode, "DELTA_TAI_GPS",  DELTA_TAI_GPS);
    writeXml(xmlNode, "DELTA_TT_GPS",   DELTA_TT_GPS);
    writeXml(xmlNode, "J2000",          STRING_J2000);

    for(UInt i=0; i<DELTA_UTC_GPS.size(); i++)
    {
      XmlNodePtr xmlNode2 = createXmlNode(xmlNode, "leapSecond");
      writeXml(xmlNode2, "MJD",            MJD_UTC_GPS.at(i));
      writeXml(xmlNode2, "DELTA_UTC_GPS",  DELTA_UTC_GPS.at(i));
    }

    OutFile file(fileName);
    XmlNode::write(file, xmlNode);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
