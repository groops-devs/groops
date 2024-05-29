/***********************************************/
/**
* @file platformSelectorEquipment.h
*
* @brief Select platforms with wildcards.
* @see PlatformSelector
*
* @author Torsten Mayer-Guerr
* @date 2024-02-24
*
*/
/***********************************************/

#ifndef __GROOPS_PLATFORMSELECTOREQUIPMENT__
#define __GROOPS_PLATFORMSELECTOREQUIPMENT__

// Latex documentation
#ifdef DOCSTRING_PlatformSelector
static const char *docstringPlatformSelectorEquipment = R"(
\subsection{Equipment}\label{platformSelectorType:equipment}
Select all platforms which has the specified equipment in the processed time interval.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "base/string.h"
#include "config/config.h"
#include "classes/platformSelector/platformSelector.h"

/***** CLASS ***********************************/

/** @brief Select platforms with wildcards.
* @ingroup platformSelectorGroup
* @see PlatformSelector */
class PlatformSelectorEquipment : public PlatformSelectorBase
{
  PlatformEquipment::Type type;
  std::regex patternName, patternSerial;
  std::regex patternRadome, patternVersion;
  std::regex patternCospar, patternNorad, patternSic, patternSp3;

public:
  PlatformSelectorEquipment(Config &config);
  void select(const Time &timeStart, const Time &timeEnd, const std::vector<const Platform*> &platforms, std::vector<Byte> &selected) const override;
};

/***********************************************/

inline PlatformSelectorEquipment::PlatformSelectorEquipment(Config &config)
{
  try
  {
    std::string choice;
    std::string name("*"), serial("*"), radome("*"), version("*");
    std::string cospar("*"), norad("*"), sic("*"), sp3("*");

    readConfig(config, "name",   name,   Config::OPTIONAL, "*", "wildcards: * and ?");
    readConfig(config, "serial", serial, Config::OPTIONAL, "*", "wildcards: * and ?");
    if(readConfigChoice(config, "equipmentType", choice, Config::MUSTSET, "", "equipment type"))
    {
      if(readConfigChoiceElement(config, "all", choice, "all types")) type = PlatformEquipment::UNDEFINED;
      if(readConfigChoiceElement(config, "gnssAntenna", choice, "antennas"))
      {
        type = PlatformEquipment::GNSSANTENNA;
        readConfig(config, "radome", radome, Config::OPTIONAL, "*", "wildcards: * and ?");
      }
      if(readConfigChoiceElement(config, "gnssReceiver", choice, "receivers"))
      {
        type = PlatformEquipment::GNSSRECEIVER;
        readConfig(config, "version", version, Config::OPTIONAL, "*", "wildcards: * and ?");
      }
      if(readConfigChoiceElement(config, "slrStation",          choice, "SLR station"))          type = PlatformEquipment::SLRSTATION;
      if(readConfigChoiceElement(config, "slrRetroReflector",   choice, "laser retroreflector")) type = PlatformEquipment::LASERRETROREFLECTOR;
      if(readConfigChoiceElement(config, "satelliteIdentifier", choice, "satellite identifier"))
      {
        type = PlatformEquipment::SATELLITEIDENTIFIER;
        readConfig(config, "cospar", cospar, Config::OPTIONAL, "*", "wildcards: * and ?");
        readConfig(config, "norad",  norad,  Config::OPTIONAL, "*", "wildcards: * and ?");
        readConfig(config, "sic",    sic,    Config::OPTIONAL, "*", "wildcards: * and ?");
        readConfig(config, "sp3",    sp3,    Config::OPTIONAL, "*", "wildcards: * and ?");
      }
      if(readConfigChoiceElement(config, "other", choice, "other types")) type = PlatformEquipment::OTHER;
      endChoice(config);
    }
    readConfig(config, "exclude", exclude, Config::DEFAULT, "0", "deselect matching platforms");
    if(isCreateSchema(config)) return;

    patternName     = String::wildcard2regex(name);
    patternSerial   = String::wildcard2regex(serial);
    patternRadome   = String::wildcard2regex(radome);
    patternVersion  = String::wildcard2regex(version);
    patternCospar   = String::wildcard2regex(cospar);
    patternNorad    = String::wildcard2regex(norad);
    patternSic      = String::wildcard2regex(sic);
    patternSp3      = String::wildcard2regex(sp3);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void PlatformSelectorEquipment::select(const Time &timeStart, const Time &timeEnd, const std::vector<const Platform*> &platforms, std::vector<Byte> &selected) const
{
  try
  {
    for(UInt i=0; i<platforms.size(); i++)
      if(platforms.at(i))
        for(const auto &eq : platforms.at(i)->equipments)
          if(((type == PlatformEquipment::UNDEFINED) || (eq->getType() == type)) &&
             (eq->timeStart < timeEnd) && (eq->timeEnd > timeStart))
          {
            if(!std::regex_match(eq->name, patternName) || !std::regex_match(eq->serial, patternSerial))
              continue;
            switch(eq->getType())
            {
              case PlatformEquipment::GNSSANTENNA:
                if(std::regex_match(std::dynamic_pointer_cast<PlatformGnssAntenna>(eq)->radome, patternRadome))
                  selected.at(i) = !exclude;
                break;
              case PlatformEquipment::GNSSRECEIVER:
                if(std::regex_match(std::dynamic_pointer_cast<PlatformGnssReceiver>(eq)->version, patternVersion))
                  selected.at(i) = !exclude;
                break;
              case PlatformEquipment::SATELLITEIDENTIFIER:
                if(std::regex_match(std::dynamic_pointer_cast<PlatformSatelliteIdentifier>(eq)->cospar, patternCospar) &&
                   std::regex_match(std::dynamic_pointer_cast<PlatformSatelliteIdentifier>(eq)->norad,  patternNorad)  &&
                   std::regex_match(std::dynamic_pointer_cast<PlatformSatelliteIdentifier>(eq)->sic,    patternSic)    &&
                   std::regex_match(std::dynamic_pointer_cast<PlatformSatelliteIdentifier>(eq)->sp3,    patternSp3))
                selected.at(i) = !exclude;
                break;
              case PlatformEquipment::SLRSTATION:          [[fallthrough]];
              case PlatformEquipment::LASERRETROREFLECTOR: [[fallthrough]];
              case PlatformEquipment::OTHER:               [[fallthrough]];
              case PlatformEquipment::UNDEFINED:           [[fallthrough]];
              default:
                selected.at(i) = !exclude;
            }
          }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
