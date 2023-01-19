/***********************************************/
/**
* @file fileGnssStationInfo.cpp
*
* @brief Description of GNSS stations. Can also be used for GNSS satellites.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2012-11-22
*
*/
/***********************************************/

#define DOCSTRING_FILEFORMAT_GnssStationInfo

#include "base/import.h"
#include "inputOutput/logging.h"
#include "inputOutput/fileArchive.h"
#include "files/fileFormatRegister.h"
#include "files/fileGnssAntennaDefinition.h"
#include "files/fileGnssStationInfo.h"

GROOPS_REGISTER_FILEFORMAT(GnssStationInfo, "gnssStationInfo")

/***********************************************/

UInt GnssStationInfo::findAntenna(const Time &time) const
{
  try
  {
    auto iter = std::lower_bound(antenna.begin(), antenna.end(), time,
                                 [](const GnssAntennaInfo &info, const Time &time) {return info.timeEnd <= time;});
    if((iter == antenna.end()) || (iter->timeStart > time))
      return NULLINDEX;
    return std::distance(antenna.begin(), iter);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

UInt GnssStationInfo::findReceiver(const Time &time) const
{
  try
  {
    auto iter = std::lower_bound(receiver.begin(), receiver.end(), time,
                                 [](const GnssReceiverInfo &info, const Time &time) {return info.timeEnd <= time;});
    if((iter == receiver.end()) || (iter->timeStart > time))
      return NULLINDEX;
    return std::distance(receiver.begin(), iter);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector GnssStationInfo::antennaVariations(const Time &time, Angle azimut, Angle elevation, const std::vector<GnssType> &types, GnssAntennaDefinition::NoPatternFoundAction noPatternFoundAction) const
{
  try
  {
    const UInt idAnt = findAntenna(time);
    if(idAnt == NULLINDEX)
      throw(Exception(markerName+"."+markerNumber+": no antenna definition found at "+time.dateTimeStr()));
    if(!antenna.at(idAnt).antennaDef)
      throw(Exception("no antenna definition for "+antenna.at(idAnt).str()));
    return antenna.at(idAnt).antennaDef->antennaVariations(azimut, elevation, types, noPatternFoundAction);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector GnssStationInfo::accuracy(const Time &time, Angle azimut, Angle elevation, const std::vector<GnssType> &types, GnssAntennaDefinition::NoPatternFoundAction noPatternFoundAction) const
{
  try
  {
    UInt idAnt = findAntenna(time);
    if(idAnt == NULLINDEX)
      throw(Exception(markerName+"."+markerNumber+": no antenna accuracy found at "+time.dateTimeStr()));
    if(!antenna.at(idAnt).accuracyDef)
      throw(Exception("no accuracy definition for "+antenna.at(idAnt).str()));
    return antenna.at(idAnt).accuracyDef->antennaVariations(azimut, elevation, types, noPatternFoundAction);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector3d GnssStationInfo::referencePoint(const Time &time) const
{
  try
  {
    if(!referencePoints.size())
      return Vector3d(0,0,0);

    auto iter = std::lower_bound(referencePoints.begin(), referencePoints.end(), time,
                                 [](const GnssReferencePointInfo &info, const Time &time) {return info.timeEnd <= time;});
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

void GnssStationInfo::fillAntennaPattern(const std::vector<GnssAntennaDefinitionPtr> &antennaList)
{
  try
  {
    for(auto &antennaInfo : antenna)
      antennaInfo.antennaDef = GnssAntennaDefinition::find(antennaList, antennaInfo.name, antennaInfo.serial, antennaInfo.radome);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssStationInfo::fillAntennaAccuracy(const std::vector<GnssAntennaDefinitionPtr> &antennaList)
{
  try
  {
    for(auto &antennaInfo : antenna)
      antennaInfo.accuracyDef = GnssAntennaDefinition::find(antennaList, antennaInfo.name, antennaInfo.serial, antennaInfo.radome);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssStationInfo::fillReceiverDefinition(const std::vector<GnssReceiverDefinitionPtr> &receiverList)
{
  try
  {
    for(auto &recv : receiver)
      recv.receiverDef = GnssReceiverDefinition::find(receiverList, recv.name, recv.serial, recv.version);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

template<> void save(OutArchive &ar, const GnssStationInfo &x)
{
  try
  {
    ar<<nameValue("markerName",     x.markerName);
    ar<<nameValue("markerNumber",   x.markerNumber);
    ar<<nameValue("comment",        x.comment);
    ar<<nameValue("approxPosition", x.approxPosition);
    ar<<nameValue("antenna",        x.antenna);
    ar<<nameValue("receiver",       x.receiver);
    ar<<nameValue("referencePoint", x.referencePoints);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void load(InArchive  &ar, GnssStationInfo &x)
{
  try
  {
    ar>>nameValue("markerName",     x.markerName);
    ar>>nameValue("markerNumber",   x.markerNumber);
    ar>>nameValue("comment",        x.comment);
    ar>>nameValue("approxPosition", x.approxPosition);
    ar>>nameValue("antenna",        x.antenna);
    ar>>nameValue("receiver",       x.receiver);
    ar>>nameValue("referencePoint", x.referencePoints);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

template<> void save(OutArchive &ar, const GnssAntennaInfo &x)
{
  try
  {
    ar<<nameValue("name",      x.name);
    ar<<nameValue("serial",    x.serial);
    ar<<nameValue("radome",    x.radome);
    ar<<nameValue("comment",   x.comment);
    ar<<nameValue("timeStart", x.timeStart);
    ar<<nameValue("timeEnd",   x.timeEnd);
    ar<<nameValue("position",  x.position);
    ar<<nameValue("local2antennaFrame",  x.local2antennaFrame);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void load(InArchive  &ar, GnssAntennaInfo &x)
{
  try
  {
    ar>>nameValue("name",      x.name);
    ar>>nameValue("serial",    x.serial);
    ar>>nameValue("radome",    x.radome);
    ar>>nameValue("comment",   x.comment);
    ar>>nameValue("timeStart", x.timeStart);
    ar>>nameValue("timeEnd",   x.timeEnd);
    ar>>nameValue("position",  x.position);
    ar>>nameValue("local2antennaFrame",  x.local2antennaFrame);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

template<> void save(OutArchive &ar, const GnssReceiverInfo &x)
{
  try
  {
    ar<<nameValue("name",      x.name);
    ar<<nameValue("serial",    x.serial);
    ar<<nameValue("version",   x.version);
    ar<<nameValue("comment",   x.comment);
    ar<<nameValue("timeStart", x.timeStart);
    ar<<nameValue("timeEnd",   x.timeEnd);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void load(InArchive  &ar, GnssReceiverInfo &x)
{
  try
  {
    ar>>nameValue("name",      x.name);
    ar>>nameValue("serial",    x.serial);
    ar>>nameValue("version",   x.version);
    ar>>nameValue("comment",   x.comment);
    ar>>nameValue("timeStart", x.timeStart);
    ar>>nameValue("timeEnd",   x.timeEnd);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

template<> void save(OutArchive &ar, const GnssReferencePointInfo &x)
{
  try
  {
    ar<<nameValue("comment",    x.comment);
    ar<<nameValue("pointStart", x.pointStart);
    ar<<nameValue("pointEnd",   x.pointEnd);
    ar<<nameValue("timeStart",  x.timeStart);
    ar<<nameValue("timeEnd",    x.timeEnd);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void load(InArchive  &ar, GnssReferencePointInfo &x)
{
  try
  {
    ar>>nameValue("comment",    x.comment);
    ar>>nameValue("pointStart", x.pointStart);
    ar>>nameValue("pointEnd",   x.pointEnd);
    ar>>nameValue("timeStart",  x.timeStart);
    ar>>nameValue("timeEnd",    x.timeEnd);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void writeFileGnssStationInfo(const FileName &fileName, const GnssStationInfo &x)
{
  try
  {
    OutFileArchive file(fileName, FILE_GNSSSTATIONINFO_TYPE);
    file<<nameValue("stationCount", 1);
    file<<nameValue("station", x);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void writeFileGnssStationInfo(const FileName &fileName, const std::vector<GnssStationInfo> &x)
{
  try
  {
    OutFileArchive file(fileName, FILE_GNSSSTATIONINFO_TYPE);
    file<<nameValue("stationCount", x.size());
    for(UInt i=0; i<x.size(); i++)
      file<<nameValue("station", x.at(i));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void readFileGnssStationInfo(const FileName &fileName, GnssStationInfo &x)
{
  try
  {
    InFileArchive file(fileName, FILE_GNSSSTATIONINFO_TYPE);
    if(file.version() < 20190304)
      throw(Exception(fileName.str()+": old GnssStationInfo file, definition of reference frames changed"));

    UInt count;
    file>>nameValue("stationCount", count);
    if(count>1)
      logWarning<<fileName<<" contain more than one station, only the first is used"<<Log::endl;
    file>>nameValue("station", x);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void readFileGnssStationInfo(const FileName &fileName, std::vector<GnssStationInfo> &x)
{
  try
  {
    InFileArchive file(fileName, FILE_GNSSSTATIONINFO_TYPE);
    if(file.version() < 20190304)
      throw(Exception(fileName.str()+": old GnssStationInfo file, definition of reference frames changed"));

    UInt count;
    file>>nameValue("stationCount", count);
    x.resize(count);
    for(UInt i=0; i<x.size(); i++)
      file>>nameValue("station", x.at(i));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/
