/***********************************************/
/**
* @file fileGnssStationInfo.h
*
* @brief Description of GNSS stations. Can also be used for GNSS satellites.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2012-11-22
*
*/
/***********************************************/

#ifndef __GROOPS_FILEGNSSSTATIONINFO__
#define __GROOPS_FILEGNSSSTATIONINFO__

// Latex documentation
#ifdef DOCSTRING_FILEFORMAT_GnssStationInfo
static const char *docstringGnssStationInfo = R"(
DEPRECATED. Use \file{Platform}{platform} instead.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "base/gnssType.h"
#include "inputOutput/fileName.h"
#include "inputOutput/fileArchive.h"

#include "files/fileGnssAntennaDefinition.h"
#include "files/fileGnssReceiverDefinition.h"

/** @addtogroup filesGroup */
/// @{

/***** CONSTANTS ********************************/

const char *const FILE_GNSSSTATIONINFO_TYPE    = "stationInfo";
constexpr UInt    FILE_GNSSSTATIONINFO_VERSION = std::max(UInt(20200123), FILE_BASE_VERSION);

/***** TYPES ***********************************/

class GnssStationInfo;
class GnssAntennaInfo;
class GnssReceiverInfo;
class GnssReferencePointInfo;
class GnssReceiverDefinition;
typedef std::shared_ptr<GnssReceiverDefinition> GnssReceiverDefinitionPtr;

/***** CLASS ***********************************/

class GnssStationInfo
{
  public:
  std::string markerName, markerNumber;
  std::string comment;
  Vector3d    approxPosition;
  std::vector<GnssAntennaInfo>        antenna;  ///< time sorted list of antennas
  std::vector<GnssReceiverInfo>       receiver; ///< time sorted list of receivers
  std::vector<GnssReferencePointInfo> referencePoints; ///< time sorted list of referencePoint

  /** @brief reference point in local frame.
  * Center of mass (CoM) for satellites. */
  Vector3d referencePoint(const Time &time) const;

  Vector antennaVariations(const Time &time, Angle azimut, Angle elevation, const std::vector<GnssType> &type, GnssAntennaDefinition::NoPatternFoundAction noPatternFoundAction) const;
  Vector accuracy         (const Time &time, Angle azimut, Angle elevation, const std::vector<GnssType> &type, GnssAntennaDefinition::NoPatternFoundAction noPatternFoundAction) const;

  void   fillAntennaPattern (const std::vector<GnssAntennaDefinitionPtr> &antennaList);
  void   fillAntennaAccuracy(const std::vector<GnssAntennaDefinitionPtr> &antennaList);
  void   fillReceiverDefinition(const std::vector<GnssReceiverDefinitionPtr> &receiverList);
  UInt   findAntenna(const Time &time) const;
  UInt   findReceiver(const Time &time) const;
};


/***** CLASS ***********************************/

class GnssAntennaInfo
{
  public:
  std::string name, serial;
  std::string radome;
  std::string comment;
  Time        timeStart, timeEnd;
  Vector3d    position;            // position of antenna reference in north, east, up or vehicle system
  Transform3d local2antennaFrame;  // north, east, up or vehicle system -> antenna system
  GnssAntennaDefinitionPtr antennaDef;
  GnssAntennaDefinitionPtr accuracyDef;

  std::string str() const {return GnssAntennaDefinition::str(name, serial, radome);}
};

/***** CLASS ***********************************/

class GnssReceiverInfo
{
  public:
  std::string name, serial;
  std::string version; // software version
  std::string comment;
  Time        timeStart, timeEnd;
  GnssReceiverDefinitionPtr receiverDef;

  std::string str() const {return GnssReceiverDefinition::str(name, serial, version);}
};

/***** CLASS ***********************************/

class GnssReferencePointInfo
{
  public:
  std::string comment;
  Vector3d    pointStart, pointEnd;
  Time        timeStart, timeEnd;
};

/***** FUNCTIONS *******************************/

template<> void save(OutArchive &ar, const GnssStationInfo &x);
template<> void save(OutArchive &ar, const GnssAntennaInfo &x);
template<> void save(OutArchive &ar, const GnssReceiverInfo &x);

template<> void load(InArchive  &ar, GnssStationInfo &x);
template<> void load(InArchive  &ar, GnssAntennaInfo &x);
template<> void load(InArchive  &ar, GnssReceiverInfo &x);

/** @brief Write into a GnssStationInfo file. */
void writeFileGnssStationInfo(const FileName &fileName, const GnssStationInfo &x);

/** @brief Write into a GnssStationInfo file. */
void writeFileGnssStationInfo(const FileName &fileName, const std::vector<GnssStationInfo> &x);

/** @brief Read from a GnssStationInfo file. */
void readFileGnssStationInfo(const FileName &fileName, GnssStationInfo &x);

/** @brief Read from a GnssStationInfo file. */
void readFileGnssStationInfo(const FileName &fileName, std::vector<GnssStationInfo> &x);

/***********************************************/

/// @}

#endif /* __GROOPS___ */
