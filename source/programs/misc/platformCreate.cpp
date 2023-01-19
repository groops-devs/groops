/***********************************************/
/**
* @file platformCreate.cpp
*
* @brief create a platform file equipped with instruments.
*
* @author Torsten Mayer-Guerr
* @date 2022-11-11
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Create a \file{Platform file}{platform} from scratch by defining attributes such as
\config{markerName}, \config{markerNumber}, \config{comment}, \config{approxPosition},
\config{equipment}.

See also \program{GnssAntex2AntennaDefinition} and \program{GnssStationLog2Platform}.

\fig{!hb}{0.8}{fileFormatPlatform}{fig:platformCreate}{Platform for stations, LEOs, and GNSS satellites.}
)";

/***********************************************/

#include "programs/program.h"
#include "files/filePlatform.h"

/***** CLASS ***********************************/

/** @brief create a platform file equipped with instruments.
* @ingroup programsGroup */
class PlatformCreate
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(PlatformCreate, SINGLEPROCESS, "create a platform file equipped with instruments", Misc)

/***********************************************/

static PlatformEquipmentPtr createEquipmentGnssAntenna(Config &config)
{
  try
  {
    auto var = new PlatformGnssAntenna();
    Angle angleX, angleY, angleZ;
    Bool  flipx, flipy, flipz;

    readConfig(config, "name",      var->name,         Config::MUSTSET,  "",  "");
    readConfig(config, "serial",    var->serial,       Config::OPTIONAL, "",  "");
    readConfig(config, "radome",    var->radome,       Config::OPTIONAL, "",  "");
    readConfig(config, "comment",   var->comment,      Config::OPTIONAL, "",  "");
    readConfig(config, "timeStart", var->timeStart,    Config::OPTIONAL, "",  "");
    readConfig(config, "timeEnd",   var->timeEnd,      Config::OPTIONAL, "",  "");
    readConfig(config, "positionX", var->position.x(), Config::MUSTSET,  "0", "[m] ARP in north, east, up or vehicle system");
    readConfig(config, "positionY", var->position.y(), Config::MUSTSET,  "0", "[m] ARP in north, east, up or vehicle system");
    readConfig(config, "positionZ", var->position.z(), Config::MUSTSET,  "0", "[m] ARP in north, east, up or vehicle system");
    readConfig(config, "rotationX", angleX,            Config::DEFAULT,  "0", "[degree] from local/vehicle to left-handed antenna system");
    readConfig(config, "rotationY", angleY,            Config::DEFAULT,  "0", "[degree] from local/vehicle to left-handed antenna system");
    readConfig(config, "rotationZ", angleZ,            Config::DEFAULT,  "0", "[degree] from local/vehicle to left-handed antenna system");
    readConfig(config, "flipX",     flipx,             Config::DEFAULT,  "0", "flip x-axis (after rotation)");
    readConfig(config, "flipY",     flipy,             Config::DEFAULT,  "0", "flip y-axis (after rotation)");
    readConfig(config, "flipZ",     flipz,             Config::DEFAULT,  "0", "flip z-axis (after rotation)");
    if(isCreateSchema(config))
      return PlatformEquipmentPtr(var);

    if(var->timeEnd == Time())
      var->timeEnd = date2time(2500,1,1);
    var->local2antennaFrame = rotaryZ(angleZ) * rotaryY(angleY) * rotaryX(angleX);
    if(flipx) var->local2antennaFrame = flipX() * var->local2antennaFrame;
    if(flipy) var->local2antennaFrame = flipY() * var->local2antennaFrame;
    if(flipz) var->local2antennaFrame = flipZ() * var->local2antennaFrame;

    return PlatformEquipmentPtr(var);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

static PlatformEquipmentPtr createEquipmentGnssReceiver(Config &config)
{
  try
  {
    auto var = new PlatformGnssReceiver();

    readConfig(config, "name",      var->name,      Config::MUSTSET,  "",  "");
    readConfig(config, "serial",    var->serial,    Config::OPTIONAL, "",  "");
    readConfig(config, "version",   var->version,   Config::OPTIONAL, "",  "");
    readConfig(config, "comment",   var->comment,   Config::OPTIONAL, "",  "");
    readConfig(config, "timeStart", var->timeStart, Config::OPTIONAL, "",  "");
    readConfig(config, "timeEnd",   var->timeEnd,   Config::OPTIONAL, "",  "");

    if(var->timeEnd == Time())
      var->timeEnd = date2time(2500,1,1);

    return PlatformEquipmentPtr(var);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

static PlatformEquipmentPtr createEquipmentOther(Config &config)
{
  try
  {
    auto var = new PlatformEquipment();

    readConfig(config, "name",      var->name,         Config::MUSTSET,  "",  "");
    readConfig(config, "serial",    var->serial,       Config::OPTIONAL, "",  "");
    readConfig(config, "comment",   var->comment,      Config::OPTIONAL, "",  "");
    readConfig(config, "timeStart", var->timeStart,    Config::OPTIONAL, "",  "");
    readConfig(config, "timeEnd",   var->timeEnd,      Config::OPTIONAL, "",  "");
    readConfig(config, "positionX", var->position.x(), Config::OPTIONAL, "0", "[m] in north, east, up or vehicle system");
    readConfig(config, "positionY", var->position.y(), Config::OPTIONAL, "0", "[m] in north, east, up or vehicle system");
    readConfig(config, "positionZ", var->position.z(), Config::OPTIONAL, "0", "[m] in north, east, up or vehicle system");

    if(var->timeEnd == Time())
      var->timeEnd = date2time(2500,1,1);

    return PlatformEquipmentPtr(var);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

static Bool readConfig(Config &config, const std::string &name, PlatformEquipmentPtr &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  try
  {
    std::string type;
    if(readConfigChoice(config, name, type, mustSet, defaultValue, annotation))
    {
      if(readConfigChoiceElement(config, "gnssAntenna",  type, ""))
        var = createEquipmentGnssAntenna(config);
      if(readConfigChoiceElement(config, "gnssReceiver", type, ""))
        var = createEquipmentGnssReceiver(config);
      if(readConfigChoiceElement(config, "other",        type, ""))
        var = createEquipmentOther(config);
      endChoice(config);
      return TRUE;
    }

    return FALSE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

static Bool readConfig(Config &config, const std::string &name, Platform::ReferencePoint &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  try
  {
    if(!readConfigSequence(config, name, mustSet, defaultValue, annotation))
      return FALSE;

    readConfig(config, "comment",   var.comment,        Config::OPTIONAL, "",  "");
    readConfig(config, "xStart",    var.pointStart.x(), Config::MUSTSET,  "0", "[m] in north, east, up or vehicle system");
    readConfig(config, "yStart",    var.pointStart.y(), Config::MUSTSET,  "0", "linear motion between start and end");
    readConfig(config, "zStart",    var.pointStart.z(), Config::MUSTSET,  "0", "");
    readConfig(config, "xEnd",      var.pointEnd.x(),   Config::MUSTSET,  "0", "[m] in north, east, up or vehicle system");
    readConfig(config, "yEnd",      var.pointEnd.y(),   Config::MUSTSET,  "0", "linear motion between start and end");
    readConfig(config, "zEnd",      var.pointEnd.z(),   Config::MUSTSET,  "0", "");
    readConfig(config, "timeStart", var.timeStart,      Config::OPTIONAL, "",  "");
    readConfig(config, "timeEnd",   var.timeEnd,        Config::OPTIONAL, "",  "");
    endSequence(config);
    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void PlatformCreate::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNamePlatform;
    Platform platform;

    readConfig(config, "outputfilePlatform", fileNamePlatform,            Config::MUSTSET,  "",  "");
    readConfig(config, "markerName",         platform.markerName,         Config::MUSTSET,  "",  "");
    readConfig(config, "markerNumber",       platform.markerNumber,       Config::OPTIONAL, "",  "");
    readConfig(config, "comment",            platform.comment,            Config::OPTIONAL, "",  "");
    readConfig(config, "approxPositionX",    platform.approxPosition.x(), Config::DEFAULT,  "0", "[m] in TRF");
    readConfig(config, "approxPositionY",    platform.approxPosition.y(), Config::DEFAULT,  "0", "[m] in TRF");
    readConfig(config, "approxPositionZ",    platform.approxPosition.z(), Config::DEFAULT,  "0", "[m] in TRF");
    readConfig(config, "equipment",          platform.equipments,         Config::OPTIONAL, "",  "");
    readConfig(config, "referencePoint",     platform.referencePoints,    Config::OPTIONAL, "",  "e.g. center of mass in satellite frame");
    if(isCreateSchema(config)) return;

    // check reference points
    // ----------------------
    if(platform.referencePoints.size())
    {
      for(UInt i=0; i<platform.referencePoints.size()-1; i++)
        if(platform.referencePoints.at(i).timeEnd == Time())
          platform.referencePoints.at(i).timeEnd = platform.referencePoints.at(i+1).timeStart;
      if(platform.referencePoints.back().timeEnd == Time())
        platform.referencePoints.back().timeEnd = date2time(2500,1,1);
    }

    logStatus<<"write platform to <"<<fileNamePlatform<<">"<<Log::endl;
    writeFilePlatform(fileNamePlatform, platform);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
