/***********************************************/
/**
* @file gnssStationInfoCreate.cpp
*
* @brief DEPRECATED since 2024-12-02. Use PlatformCreate instead.
*
* @author Torsten Mayer-Guerr
* @date 2012-11-24
*
* @deprecated Use PlatformCreate instead.
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
DEPRECATED since 2024-12-02. Please use \program{PlatformCreate} instead.
)";

/***********************************************/

#include "programs/program.h"
#include "files/filePlatform.h"

/***** CLASS ***********************************/

/** @brief DEPRECATED since 2024-12-02. Please use PlatformCreate instead.
* @ingroup programsGroup */
class GnssStationInfoCreate
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GnssStationInfoCreate, SINGLEPROCESS, "DEPRECATED since 2024-12-02. Please use PlatformCreate instead.", Deprecated)
GROOPS_RENAMED_PROGRAM(GnssCreateStationInfo, GnssStationInfoCreate, date2time(2019, 9, 5))

/***********************************************/

static Bool readConfig(Config &config, const std::string &name, std::shared_ptr<PlatformGnssAntenna> &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  try
  {
    if(!readConfigSequence(config, name, mustSet, defaultValue, annotation))
      return FALSE;

    var = std::make_shared<PlatformGnssAntenna>();
    Angle angleX, angleY, angleZ;
    Bool  flipx, flipy, flipz;

    readConfig(config, "name",        var->name,         Config::MUSTSET,  "",  "");
    readConfig(config, "serial",      var->serial,       Config::OPTIONAL, "",  "");
    readConfig(config, "radome",      var->radome,       Config::OPTIONAL, "",  "");
    readConfig(config, "comment",     var->comment,      Config::OPTIONAL, "",  "");
    readConfig(config, "timeStart",   var->timeStart,    Config::OPTIONAL, "",  "");
    readConfig(config, "timeEnd",     var->timeEnd,      Config::OPTIONAL, "",  "");
    readConfig(config, "positionX",   var->position.x(), Config::MUSTSET,  "0", "[m] ARP in north, east, up or vehicle system");
    readConfig(config, "positionY",   var->position.y(), Config::MUSTSET,  "0", "[m] ARP in north, east, up or vehicle system");
    readConfig(config, "positionZ",   var->position.z(), Config::MUSTSET,  "0", "[m] ARP in north, east, up or vehicle system");
    readConfig(config, "rotationX",   angleX,            Config::DEFAULT,  "0", "[degree] from local/vehicle to left-handed antenna system");
    readConfig(config, "rotationY",   angleY,            Config::DEFAULT,  "0", "[degree] from local/vehicle to left-handed antenna system");
    readConfig(config, "rotationZ",   angleZ,            Config::DEFAULT,  "0", "[degree] from local/vehicle to left-handed antenna system");
    readConfig(config, "flipX",       flipx,             Config::DEFAULT,  "0", "flip x-axis (after rotation)");
    readConfig(config, "flipY",       flipy,             Config::DEFAULT,  "0", "flip y-axis (after rotation)");
    readConfig(config, "flipZ",       flipz,             Config::DEFAULT,  "0", "flip z-axis (after rotation)");
    endSequence(config);
    if(isCreateSchema(config))
      return TRUE;

    if(var->timeEnd == Time())
      var->timeEnd = date2time(2500,1,1);
    var->local2antennaFrame = rotaryZ(angleZ) * rotaryY(angleY) * rotaryX(angleX);
    if(flipx) var->local2antennaFrame = flipX() * var->local2antennaFrame;
    if(flipy) var->local2antennaFrame = flipY() * var->local2antennaFrame;
    if(flipz) var->local2antennaFrame = flipZ() * var->local2antennaFrame;

    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

static Bool readConfig(Config &config, const std::string &name, std::shared_ptr<PlatformGnssReceiver> &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  try
  {
    if(!readConfigSequence(config, name, mustSet, defaultValue, annotation))
      return FALSE;

    var = std::make_shared<PlatformGnssReceiver>();

    readConfig(config, "name",      var->name,      Config::MUSTSET,  "",  "");
    readConfig(config, "serial",    var->serial,    Config::OPTIONAL, "",  "");
    readConfig(config, "version",   var->version,   Config::OPTIONAL, "",  "");
    readConfig(config, "comment",   var->comment,   Config::OPTIONAL, "",  "");
    readConfig(config, "timeStart", var->timeStart, Config::OPTIONAL, "",  "");
    readConfig(config, "timeEnd",   var->timeEnd,   Config::OPTIONAL, "",  "");
    endSequence(config);
    return TRUE;
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

void GnssStationInfoCreate::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName outFileName;
    Platform station;
    std::vector<std::shared_ptr<PlatformGnssAntenna>>  antennas;
    std::vector<std::shared_ptr<PlatformGnssReceiver>> receivers;

    readConfig(config, "outputfileStationInfo", outFileName,                Config::MUSTSET,  "",  "");
    readConfig(config, "markerName",            station.markerName,         Config::MUSTSET,  "",  "");
    readConfig(config, "markerNumber",          station.markerNumber,       Config::OPTIONAL, "",  "");
    readConfig(config, "comment",               station.comment,            Config::OPTIONAL, "",  "");
    readConfig(config, "approxPositionX",       station.approxPosition.x(), Config::DEFAULT,  "0", "[m] in TRF");
    readConfig(config, "approxPositionY",       station.approxPosition.y(), Config::DEFAULT,  "0", "[m] in TRF");
    readConfig(config, "approxPositionZ",       station.approxPosition.z(), Config::DEFAULT,  "0", "[m] in TRF");
    readConfig(config, "antenna",               antennas,                   Config::MUSTSET,  "",  "");
    readConfig(config, "receiver",              receivers,                  Config::OPTIONAL, "",  "");
    readConfig(config, "referencePoint",        station.referencePoints,    Config::OPTIONAL, "",  "e.g. center of mass in satellite frame");
    if(isCreateSchema(config)) return;

    logWarning<<"DEPRECATED since 2024-12-02. Please use PlatformCreate instead."<<Log::endl;

    // ============================

    // check antennas
    // --------------
    if(antennas.size())
    {
      for(UInt i=0; i<antennas.size()-1; i++)
        if(antennas.at(i)->timeEnd == Time())
          antennas.at(i)->timeEnd = antennas.at(i+1)->timeStart;
      if(antennas.back()->timeEnd == Time())
        antennas.back()->timeEnd = date2time(2500,1,1);
      station.equipments.insert(station.equipments.end(), antennas.begin(), antennas.end());
    }

    // check receivers
    // ----------------
    if(receivers.size())
    {
      for(UInt i=0; i<receivers.size()-1; i++)
        if(receivers.at(i)->timeEnd == Time())
          receivers.at(i)->timeEnd = receivers.at(i+1)->timeStart;
      if(receivers.back()->timeEnd == Time())
        receivers.back()->timeEnd = date2time(2500,1,1);
      station.equipments.insert(station.equipments.end(), receivers.begin(), receivers.end());
    }

    // check reference points
    // ----------------------
    if(station.referencePoints.size())
    {
      for(UInt i=0; i<station.referencePoints.size()-1; i++)
        if(station.referencePoints.at(i).timeEnd == Time())
          station.referencePoints.at(i).timeEnd = station.referencePoints.at(i+1).timeStart;
      if(station.referencePoints.back().timeEnd == Time())
        station.referencePoints.back().timeEnd = date2time(2500,1,1);
    }

    // ============================

    logStatus<<"write platform <"<<outFileName<<">"<<Log::endl;
    writeFilePlatform(outFileName, station);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/
