/***********************************************/
/**
* @file gnssStationInfoCreate.cpp
*
* @brief create station information file for GNSS receivers.
*
* @author Torsten Mayer-Guerr
* @date 2012-11-24
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Create a \file{GnssStationInfo file}{gnssStationInfo} from scratch by defining attributes such as
\config{markerName}, \config{markerNumber}, \config{comment}, \config{approxPosition},
\config{antenna} and \config{receiver}.

See also \program{GnssAntex2AntennaDefinition} and \program{GnssStationLog2StationInfo}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileGnssStationInfo.h"

/***** CLASS ***********************************/

/** @brief create station information file for GNSS receivers.
* @ingroup programsGroup */
class GnssStationInfoCreate
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GnssStationInfoCreate, SINGLEPROCESS, "create station information file for GNSS receivers", Gnss)
GROOPS_RENAMED_PROGRAM(GnssCreateStationInfo, GnssStationInfoCreate, date2time(2019, 9, 5))

/***********************************************/

static Bool readConfig(Config &config, const std::string &name, GnssAntennaInfo &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  try
  {
    if(!readConfigSequence(config, name, mustSet, defaultValue, annotation))
      return FALSE;

    Angle angleX, angleY, angleZ;
    Bool  flipx, flipy, flipz;

    readConfig(config, "name",        var.name,         Config::MUSTSET,  "",  "");
    readConfig(config, "serial",      var.serial,       Config::OPTIONAL, "",  "");
    readConfig(config, "radome",      var.radome,       Config::OPTIONAL, "",  "");
    readConfig(config, "comment",     var.comment,      Config::OPTIONAL, "",  "");
    readConfig(config, "timeStart",   var.timeStart,    Config::OPTIONAL, "",  "");
    readConfig(config, "timeEnd",     var.timeEnd,      Config::OPTIONAL, "",  "");
    readConfig(config, "positionX",   var.position.x(), Config::MUSTSET,  "0", "[m] ARP in north, east, up or vehicle system");
    readConfig(config, "positionY",   var.position.y(), Config::MUSTSET,  "0", "[m] ARP in north, east, up or vehicle system");
    readConfig(config, "positionZ",   var.position.z(), Config::MUSTSET,  "0", "[m] ARP in north, east, up or vehicle system");
    readConfig(config, "rotationX",   angleX,           Config::DEFAULT,  "0", "[degree] from local/vehicle to left-handed antenna system");
    readConfig(config, "rotationY",   angleY,           Config::DEFAULT,  "0", "[degree] from local/vehicle to left-handed antenna system");
    readConfig(config, "rotationZ",   angleZ,           Config::DEFAULT,  "0", "[degree] from local/vehicle to left-handed antenna system");
    readConfig(config, "flipX",       flipx,            Config::DEFAULT,  "0", "flip x-axis (after rotation)");
    readConfig(config, "flipY",       flipy,            Config::DEFAULT,  "0", "flip y-axis (after rotation)");
    readConfig(config, "flipZ",       flipz,            Config::DEFAULT,  "0", "flip z-axis (after rotation)");
    endSequence(config);
    if(isCreateSchema(config))
      return TRUE;

    if(var.timeEnd == Time())
      var.timeEnd = date2time(2500,1,1);
    var.local2antennaFrame = rotaryZ(angleZ) * rotaryY(angleY) * rotaryX(angleX);
    if(flipx) var.local2antennaFrame = flipX() * var.local2antennaFrame;
    if(flipy) var.local2antennaFrame = flipY() * var.local2antennaFrame;
    if(flipz) var.local2antennaFrame = flipZ() * var.local2antennaFrame;

    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

static Bool readConfig(Config &config, const std::string &name, GnssReceiverInfo &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  try
  {
    if(!readConfigSequence(config, name, mustSet, defaultValue, annotation))
      return FALSE;

    readConfig(config, "name",      var.name,      Config::MUSTSET,  "",  "");
    readConfig(config, "serial",    var.serial,    Config::OPTIONAL, "",  "");
    readConfig(config, "version",   var.version,   Config::OPTIONAL, "",  "");
    readConfig(config, "comment",   var.comment,   Config::OPTIONAL, "",  "");
    readConfig(config, "timeStart", var.timeStart, Config::OPTIONAL, "",  "");
    readConfig(config, "timeEnd",   var.timeEnd,   Config::OPTIONAL, "",  "");
    endSequence(config);
    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

static Bool readConfig(Config &config, const std::string &name, GnssReferencePointInfo &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
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
    FileName        outFileName;
    GnssStationInfo station;

    readConfig(config, "outputfileStationInfo", outFileName,                Config::MUSTSET,  "",  "");
    readConfig(config, "markerName",            station.markerName,         Config::MUSTSET,  "",  "");
    readConfig(config, "markerNumber",          station.markerNumber,       Config::OPTIONAL, "",  "");
    readConfig(config, "comment",               station.comment,            Config::OPTIONAL, "",  "");
    readConfig(config, "approxPositionX",       station.approxPosition.x(), Config::DEFAULT,  "0", "[m] in TRF");
    readConfig(config, "approxPositionY",       station.approxPosition.y(), Config::DEFAULT,  "0", "[m] in TRF");
    readConfig(config, "approxPositionZ",       station.approxPosition.z(), Config::DEFAULT,  "0", "[m] in TRF");
    readConfig(config, "antenna",               station.antenna,            Config::MUSTSET,  "",  "");
    readConfig(config, "receiver",              station.receiver,           Config::OPTIONAL, "",  "");
    readConfig(config, "referencePoint",        station.referencePoints,    Config::OPTIONAL, "",  "e.g. center of mass in satellite frame");
    if(isCreateSchema(config)) return;

    // ============================

    // check antennas
    // --------------
    if(station.antenna.size())
    {
      for(UInt i=0; i<station.antenna.size()-1; i++)
        if(station.antenna.at(i).timeEnd == Time())
          station.antenna.at(i).timeEnd = station.antenna.at(i+1).timeStart;
      if(station.antenna.back().timeEnd == Time())
        station.antenna.back().timeEnd = date2time(2500,1,1);
    }

    // check receivers
    // ----------------
    if(station.receiver.size())
    {
      for(UInt i=0; i<station.receiver.size()-1; i++)
        if(station.receiver.at(i).timeEnd == Time())
          station.receiver.at(i).timeEnd = station.receiver.at(i+1).timeStart;
      if(station.receiver.back().timeEnd == Time())
        station.receiver.back().timeEnd = date2time(2500,1,1);
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

    logStatus<<"write antenna definition <"<<outFileName<<">"<<Log::endl;
    writeFileGnssStationInfo(outFileName, station);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/
