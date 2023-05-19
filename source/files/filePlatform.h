/***********************************************/
/**
* @file filePlatform.h
*
* @brief Platform equipped with instruments.
*
* @author Torsten Mayer-Guerr
* @date 2022-11-07
*
*/
/***********************************************/

#ifndef __GROOPS_FILEPLATFORM__
#define __GROOPS_FILEPLATFORM__

// Latex documentation
#ifdef DOCSTRING_FILEFORMAT_Platform
static const char *docstringPlatform = R"(
Defines a platform with a local coordinate frame equipped with instruments.
The platform might be a reference station, a low Earth satellite,
or a transmitting GNSS satellite and is referenced by a marker name and number.
The reference point (marker or center of mass (CoM)) can change in time
relative to the local frame.

Each equipped instrument is described at least by the following information
\begin{itemize}
\item name
\item serial number
\item coordinates in the local frame
\item a time interval in which the instrument was active
\item the orientation for antennas and reflectors.
\end{itemize}

For GNSS satellites the platform defines the PRN. The different assigned SVNs
are defined by the transmitting antennas.

Platforms for GNSS stations can be created from station log files with
\program{GnssStationLog2Platform}. Platforms for GNSS satellites
can be created from an ANTEX file with \program{GnssAntex2AntennaDefinition}.

See also \program{PlatformCreate}.

\fig{!hb}{0.8}{fileFormatPlatform}{fig:fileFormatPlatform}{Platform for stations, LEOs, and GNSS satellites.}
)";
#endif

/***********************************************/

#include "base/import.h"
#include "inputOutput/fileName.h"
#include "inputOutput/archive.h"
#include "files/fileGnssAntennaDefinition.h"
#include "files/fileGnssReceiverDefinition.h"

/** @addtogroup filesGroup */
/// @{

/***** CONSTANTS ********************************/

const char *const FILE_PLATFORM_TYPE = "platform";

/***** TYPES ***********************************/

class Platform;
class PlatformEquipment;
typedef std::shared_ptr<PlatformEquipment> PlatformEquipmentPtr;

class GnssAntennaDefinition;
class GnssReceiverDefinition;
typedef std::shared_ptr<GnssAntennaDefinition>  GnssAntennaDefinitionPtr;
typedef std::shared_ptr<GnssReceiverDefinition> GnssReceiverDefinitionPtr;

/***** CLASS ***********************************/

class Platform
{
public:
  class ReferencePoint
  {
  public:
    std::string comment;
    Vector3d    pointStart, pointEnd;
    Time        timeStart, timeEnd;
  };

  std::string                       name;
  std::string                       markerName, markerNumber;
  std::string                       comment;
  Vector3d                          approxPosition;
  std::vector<ReferencePoint>       referencePoints; ///< time sorted list of referencePoint
  std::vector<PlatformEquipmentPtr> equipments;

  /** @brief reference point in local frame.
  * Center of mass (CoM) for satellites. */
  Vector3d referencePoint(const Time &time) const;

  template<typename T> std::shared_ptr<T> findEquipment(const Time &time) const;

  void fillGnssAntennaDefinition (const std::vector<GnssAntennaDefinitionPtr> &antennaList);
  void fillGnssAccuracyDefinition(const std::vector<GnssAntennaDefinitionPtr> &antennaList);
  void fillGnssReceiverDefinition(const std::vector<GnssReceiverDefinitionPtr> &receiverList);
};

/***** CLASS ***********************************/

class PlatformEquipment
{
public:
  enum Type : Int {UNDEFINED    = 0,
                   OTHER        = 1,
                   GNSSANTENNA  = 2,
                   GNSSRECEIVER = 3,
                   SLRSTATION   = 4};

  static constexpr Type TYPE = OTHER;
  std::string comment;
  std::string name, serial;
  Time        timeStart, timeEnd;
  Vector3d    position;   // position of instrument in north, east, up or vehicle system

  /** @brief Create Equipment of given type (as shared_ptr). */
  static PlatformEquipmentPtr create(Type type);

  /** @brief Data type (e.g. GNSSANTENNA, SLRREFLECTOR). */
  virtual Type getType() const {return TYPE;}

  virtual std::string str() const {return name+"|"+serial;};

  virtual void save(OutArchive &oa) const;
  virtual void load(InArchive  &ia);
};

/***** CLASS ***********************************/

class PlatformGnssAntenna : public PlatformEquipment
{
public:
  static constexpr Type TYPE = GNSSANTENNA;
  std::string              radome;
  Transform3d              local2antennaFrame;  // north, east, up or vehicle system -> antenna system
  GnssAntennaDefinitionPtr antennaDef;
  GnssAntennaDefinitionPtr accuracyDef;

  Type getType() const override {return TYPE;}
  std::string str() const override {return GnssAntennaDefinition::str(name, serial, radome);}
  void save(OutArchive &oa) const override;
  void load(InArchive  &ia) override;
};

/***** CLASS ***********************************/

class PlatformGnssReceiver : public PlatformEquipment
{
public:
  static constexpr Type TYPE = GNSSRECEIVER;
  std::string               version; // software version
  GnssReceiverDefinitionPtr receiverDef;

  Type getType() const override {return TYPE;}
  std::string str() const override {return GnssReceiverDefinition::str(name, serial, version);}
  void save(OutArchive &oa) const override;
  void load(InArchive  &ia) override;
};

/***** CLASS ***********************************/

class PlatformSlrStation : public PlatformEquipment
{
public:
  static constexpr Type TYPE = SLRSTATION;
  Type getType() const override {return TYPE;}
  std::string str() const override {return name;}
  void save(OutArchive &oa) const override;
  void load(InArchive  &ia) override;
};

/***** FUNCTIONS *******************************/

template<> void save(OutArchive &ar, const Platform &x);
template<> void load(InArchive  &ar, Platform &x);

/** @brief Write into a Platform file. */
void writeFilePlatform(const FileName &fileName, const Platform &x);

/** @brief Read from a Platform file. */
void readFilePlatform(const FileName &fileName, Platform &x);

/***********************************************/
/***** INLINES   *******************************/
/***********************************************/

template<typename T> inline std::shared_ptr<T> Platform::findEquipment(const Time &time) const
{
  try
  {
    auto iter = std::find_if(equipments.begin(), equipments.end(), [&](const auto &x)
                            {return (x->getType() == T::TYPE) && (x->timeStart <= time) && (time < x->timeEnd);});
    if(iter == equipments.end())
      return std::shared_ptr<T>(nullptr);
    return std::dynamic_pointer_cast<T>(*iter);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

/// @}

#endif /* __GROOPS___ */
