/***********************************************/
/**
* @file filePlatform.cpp
*
* @brief Platform equipped with instruments.
*
* @author Torsten Mayer-Guerr
* @date 2022-11-07
*
*/
/***********************************************/

#define DOCSTRING_FILEFORMAT_Platform

#include "base/import.h"
#include "inputOutput/logging.h"
#include "inputOutput/fileArchive.h"
#include "files/fileFormatRegister.h"
#include "files/fileGnssAntennaDefinition.h"
#include "files/fileGnssReceiverDefinition.h"
#include "files/filePlatform.h"

GROOPS_REGISTER_FILEFORMAT(Platform, "platform")

/***********************************************/

Vector3d Platform::referencePoint(const Time &time) const
{
  try
  {
    if(!referencePoints.size())
      return Vector3d(0,0,0);

    auto iter = std::lower_bound(referencePoints.begin(), referencePoints.end(), time,
                                 [](const auto &x, const Time &time) {return x.timeEnd <= time;});
    if((iter == referencePoints.end()) || (iter->timeStart > time))
      throw(Exception(markerName+"."+markerNumber+": no reference point found at "+time.dateTimeStr()));

    // linear interpolation
    const Double tau = (time-iter->timeStart).mjd()/(iter->timeEnd-iter->timeStart).mjd();
    return (1-tau) * iter->pointStart + tau * iter->pointEnd;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Platform::fillGnssAntennaDefinition(const std::vector<GnssAntennaDefinitionPtr> &antennaList)
{
  try
  {
    for(const auto &eq : equipments)
    {
      auto antenna = std::dynamic_pointer_cast<PlatformGnssAntenna>(eq);
      if(antenna)
        antenna->antennaDef = GnssAntennaDefinition::find(antennaList, antenna->name, antenna->serial, antenna->radome);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Platform::fillGnssAccuracyDefinition(const std::vector<GnssAntennaDefinitionPtr> &antennaList)
{
  try
  {
    for(const auto &eq : equipments)
    {
      auto antenna = std::dynamic_pointer_cast<PlatformGnssAntenna>(eq);
      if(antenna)
        antenna->accuracyDef = GnssAntennaDefinition::find(antennaList, antenna->name, antenna->serial, antenna->radome);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Platform::fillGnssReceiverDefinition(const std::vector<GnssReceiverDefinitionPtr> &receiverList)
{
  try
  {
    for(const auto &eq : equipments)
    {
      auto receiver = std::dynamic_pointer_cast<PlatformGnssReceiver>(eq);
      if(receiver)
        receiver->receiverDef = GnssReceiverDefinition::find(receiverList, receiver->name, receiver->serial, receiver->version);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

PlatformEquipmentPtr PlatformEquipment::create(Type type)
{
  try
  {
    switch(type)
    {
      case OTHER:               return std::make_shared<PlatformEquipment>();
      case GNSSANTENNA:         return std::make_shared<PlatformGnssAntenna>();
      case GNSSRECEIVER:        return std::make_shared<PlatformGnssReceiver>();
      case SLRSTATION:          return std::make_shared<PlatformSlrStation>();
      case LASERRETROREFLECTOR: return std::make_shared<PlatformLaserRetroReflector>();
      case SATELLITEIDENTIFIER: return std::make_shared<PlatformSatelliteIdentifier>();
      case UNDEFINED:           break;
    }

    throw(Exception("unknown equipment type ("+static_cast<Int>(type)%"%i)"s));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void PlatformEquipment::save(OutArchive &/*ar*/) const {}
void PlatformEquipment::load(InArchive &/*ar*/) {}

/***********************************************/
/***********************************************/

void PlatformGnssAntenna::save(OutArchive &ar) const
{
  ar<<nameValue("radome",             radome);
  ar<<nameValue("local2antennaFrame", local2antennaFrame);
}

/***********************************************/

void PlatformGnssAntenna::load(InArchive &ar)
{
  ar>>nameValue("radome",             radome);
  ar>>nameValue("local2antennaFrame", local2antennaFrame);
}

/***********************************************/
/***********************************************/

void PlatformGnssReceiver::save(OutArchive &ar) const
{
  ar<<nameValue("version", version);
}

/***********************************************/

void PlatformGnssReceiver::load(InArchive &ar)
{
  ar>>nameValue("version", version);
}

/***********************************************/
/***********************************************/

void PlatformSlrStation::save(OutArchive &/*ar*/) const
{
}

/***********************************************/

void PlatformSlrStation::load(InArchive &/*ar*/)
{
}

/***********************************************/
/***********************************************/

void PlatformLaserRetroReflector::save(OutArchive &ar) const
{
  ar<<nameValue("platform2reflectorFrame", platform2reflectorFrame);
  ar<<nameValue("dZenit",                  dZenit);
  ar<<nameValue("range",                   range);
}

/***********************************************/

void PlatformLaserRetroReflector::load(InArchive &ar)
{
  ar>>nameValue("platform2reflectorFrame", platform2reflectorFrame);
  ar>>nameValue("dZenit",                  dZenit);
  ar>>nameValue("range",                   range);
}

/***********************************************/
/***********************************************/

void PlatformSatelliteIdentifier::save(OutArchive &ar) const
{
  ar<<nameValue("cospar", cospar);
  ar<<nameValue("norad",  norad);
  ar<<nameValue("sic",    sic);
  ar<<nameValue("sp3",    sp3);
}

/***********************************************/

void PlatformSatelliteIdentifier::load(InArchive &ar)
{
  ar>>nameValue("cospar", cospar);
  ar>>nameValue("norad",  norad);
  ar>>nameValue("sic",    sic);
  ar>>nameValue("sp3",    sp3);
}

/***********************************************/
/***********************************************/

template<> void save(OutArchive &ar, const Platform &x)
{
  try
  {
    ar<<nameValue("markerName",     x.markerName);
    ar<<nameValue("markerNumber",   x.markerNumber);
    ar<<nameValue("comment",        x.comment);
    ar<<nameValue("approxPosition", x.approxPosition);

    ar<<nameValue("equipmentCount", x.equipments.size());
    for(UInt i=0; i<x.equipments.size(); i++)
    {
      ar<<beginGroup("equipment");
      ar<<nameValue("type", static_cast<Int>(x.equipments.at(i)->getType()));
      ar<<nameValue("comment",   x.equipments.at(i)->comment);
      ar<<nameValue("name",      x.equipments.at(i)->name);
      ar<<nameValue("serial",    x.equipments.at(i)->serial);
      ar<<nameValue("timeStart", x.equipments.at(i)->timeStart);
      ar<<nameValue("timeEnd",   x.equipments.at(i)->timeEnd);
      ar<<nameValue("position",  x.equipments.at(i)->position);
      x.equipments.at(i)->save(ar);
      ar<<endGroup("equipment");
    }

    ar<<nameValue("referencePointCount", x.referencePoints.size());
    for(UInt i=0; i<x.referencePoints.size(); i++)
    {
      ar<<beginGroup("referencePoint");
      ar<<nameValue("comment",    x.referencePoints.at(i).comment);
      ar<<nameValue("timeStart",  x.referencePoints.at(i).timeStart);
      ar<<nameValue("timeEnd",    x.referencePoints.at(i).timeEnd);
      ar<<nameValue("pointStart", x.referencePoints.at(i).pointStart);
      ar<<nameValue("pointEnd",   x.referencePoints.at(i).pointEnd);
      ar<<endGroup("referencePoint");
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void load(InArchive &ar, Platform &x)
{
  try
  {
    ar>>nameValue("markerName",     x.markerName);
    ar>>nameValue("markerNumber",   x.markerNumber);
    ar>>nameValue("comment",        x.comment);
    ar>>nameValue("approxPosition", x.approxPosition);

    UInt equipmentCount;
    ar>>nameValue("equipmentCount", equipmentCount);
    x.equipments.resize(equipmentCount);
    for(UInt i=0; i<equipmentCount; i++)
    {
      Int typeInt;
      ar>>beginGroup("equipment");
      ar>>nameValue("type", typeInt);
      x.equipments.at(i) = PlatformEquipment::create(static_cast<PlatformEquipment::Type>(typeInt));
      ar>>nameValue("comment",   x.equipments.at(i)->comment);
      ar>>nameValue("name",      x.equipments.at(i)->name);
      ar>>nameValue("serial",    x.equipments.at(i)->serial);
      ar>>nameValue("timeStart", x.equipments.at(i)->timeStart);
      ar>>nameValue("timeEnd",   x.equipments.at(i)->timeEnd);
      ar>>nameValue("position",  x.equipments.at(i)->position);
      x.equipments.at(i)->load(ar);
      ar>>endGroup("equipment");
    }

    UInt referencePointCount;
    ar>>nameValue("referencePointCount", referencePointCount);
    x.referencePoints.resize(referencePointCount);
    for(UInt i=0; i<referencePointCount; i++)
    {
      ar>>beginGroup("referencePoint");
      ar>>nameValue("comment",    x.referencePoints.at(i).comment);
      ar>>nameValue("timeStart",  x.referencePoints.at(i).timeStart);
      ar>>nameValue("timeEnd",    x.referencePoints.at(i).timeEnd);
      ar>>nameValue("pointStart", x.referencePoints.at(i).pointStart);
      ar>>nameValue("pointEnd",   x.referencePoints.at(i).pointEnd);
      ar>>endGroup("referencePoint");
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void writeFilePlatform(const FileName &fileName, const Platform &x)
{
  try
  {
    OutFileArchive file(fileName, FILE_PLATFORM_TYPE, FILE_PLATFORM_VERSION);
    file<<nameValue("platform", x);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void readFilePlatform(const FileName &fileName, Platform &x)
{
  try
  {
    InFileArchive file(fileName, ""/*arbitrary type*/, FILE_PLATFORM_VERSION);
    if(file.type() == FILE_PLATFORM_TYPE)
      file>>nameValue("platform", x);
    else if(file.type() == "stationInfo") // old deprecated file format
    {
      UInt count;
      file>>nameValue("stationCount", count);
      if(count>1)
        logWarning<<fileName<<" contain more than one station, only the first is used"<<Log::endl;
      file>>beginGroup("station");
      file>>nameValue("markerName",     x.markerName);
      file>>nameValue("markerNumber",   x.markerNumber);
      file>>nameValue("comment",        x.comment);
      file>>nameValue("approxPosition", x.approxPosition);
      file>>beginGroup("antenna");
      {
        UInt count;
        file>>nameValue("count", count);
        for(UInt i=0; i<count; i++)
        {
          auto var = std::make_shared<PlatformGnssAntenna>();
          file>>beginGroup("cell");
          file>>nameValue("name",               var->name);
          file>>nameValue("serial",             var->serial);
          file>>nameValue("radome",             var->radome);
          file>>nameValue("comment",            var->comment);
          file>>nameValue("timeStart",          var->timeStart);
          file>>nameValue("timeEnd",            var->timeEnd);
          file>>nameValue("position",           var->position);
          file>>nameValue("local2antennaFrame", var->local2antennaFrame);
          file>>endGroup("cell");
          x.equipments.push_back(var);
        }
      }
      file>>endGroup("antenna");
      file>>beginGroup("receiver");
      {
        UInt count;
        file>>nameValue("count", count);
        for(UInt i=0; i<count; i++)
        {
          auto var = std::make_shared<PlatformGnssReceiver>();
          file>>beginGroup("cell");
          file>>nameValue("name",      var->name);
          file>>nameValue("serial",    var->serial);
          file>>nameValue("version",   var->version);
          file>>nameValue("comment",   var->comment);
          file>>nameValue("timeStart", var->timeStart);
          file>>nameValue("timeEnd",   var->timeEnd);
          file>>endGroup("cell");
          x.equipments.push_back(var);
        }
      }
      file>>endGroup("receiver");
      file>>beginGroup("referencePoint");
      {
        UInt count;
        file>>nameValue("count", count);
        x.referencePoints.resize(count);
        for(UInt i=0; i<count; i++)
        {
          file>>beginGroup("cell");
          file>>nameValue("comment",    x.referencePoints.at(i).comment);
          file>>nameValue("pointStart", x.referencePoints.at(i).pointStart);
          file>>nameValue("pointEnd",   x.referencePoints.at(i).pointEnd);
          file>>nameValue("timeStart",  x.referencePoints.at(i).timeStart);
          file>>nameValue("timeEnd",    x.referencePoints.at(i).timeEnd);
          file>>endGroup("cell");
        }
      }
      file>>endGroup("referencePoint");
      file>>endGroup("station");
    }
    else
      throw(Exception("file type is '"+file.type()+"' but must be '"+FILE_PLATFORM_TYPE+"' or 'stationInfo'"));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
