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
Information file for a station or satellite equipped with a GNSS receiver
or for a GNSS satellite.

For stations and receiving satellites it defines the time intervals and
local coordinates of installed antennas and receivers.
For GNSS satellites it defines the PRN-SVN mapping and antenna offsets
and orientations.

GNSS station infos can be created from station log files with
\program{GnssStationLog2StationInfo}. GNSS satellite transmitter infos
can be created from an ANTEX file with \program{GnssAntex2AntennaDefinition}.

See \program{GnssStationInfoCreate}.

\begin{verbatim}
<?xml version="1.0" encoding="UTF-8"?>
<groops type="stationInfo" version="20190429">
    <stationCount>1</stationCount>
    <station>
        <markerName>GRAZ</markerName>
        <markerNumber>11001M002</markerNumber>
        <comment>Graz, AUSTRIA</comment>
        <approxPosition>
            <x>4.19442412400000e+06</x>
            <y>1.16270246000000e+06</y>
            <z>4.64724519300000e+06</z>
        </approxPosition>
        <antenna>
            <count>8</count>
            ...
            <cell>
                <name>LEIAR25.R4</name>
                <serial>726444</serial>
                <radome>LEIT</radome>
                <comment/>
                <timeStart>58080.4041666666666650</timeStart>
                <timeEnd>234166.0000000000000000</timeEnd>
                <position>
                    <x>0.00000000000000e+00</x>
                    <y>0.00000000000000e+00</y>
                    <z>1.96400000000000e+00</z>
                </position>
                <local2antennaFrame>
                    <xx>1.00000000000000e+00</xx>
                    <xy>0.00000000000000e+00</xy>
                    <xz>0.00000000000000e+00</xz>
                    <yx>0.00000000000000e+00</yx>
                    <yy>1.00000000000000e+00</yy>
                    <yz>0.00000000000000e+00</yz>
                    <zx>0.00000000000000e+00</zx>
                    <zy>0.00000000000000e+00</zy>
                    <zz>1.00000000000000e+00</zz>
                </local2antennaFrame>
            </cell>
        </antenna>
        <receiver>
            <count>25</count>
            ...
            <cell>
                <name>SEPT POLARX5</name>
                <serial>4501501</serial>
                <version>5.2.0</version>
                <comment/>
                <timeStart>58903.2958333333333343</timeStart>
                <timeEnd>234166.0000000000000000</timeEnd>
            </cell>
        </receiver>
        <referencePoint>
            <count>0</count>
        </referencePoint>
    </station>
</groops>
\end{verbatim}
)";
#endif

/***********************************************/

#include "base/import.h"
#include "base/gnssType.h"
#include "inputOutput/fileName.h"
#include "inputOutput/archive.h"
#include "files/fileGnssAntennaDefinition.h"
#include "files/fileGnssReceiverDefinition.h"

/** @addtogroup filesGroup */
/// @{

/***** CONSTANTS ********************************/

const char *const FILE_GNSSSTATIONINFO_TYPE = "stationInfo";

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
