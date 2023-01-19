/***********************************************/
/**
* @file loopFileGnssStationInfo.h
*
* @brief DEPRECATDED. Use LoopPlatformEquipment instead.
*
* @author Sebastian Strasser
* @date 2018-01-29
*
*/
/***********************************************/

#ifndef __GROOPS_LOOPFILEGNSSSTATIONINFO__
#define __GROOPS_LOOPFILEGNSSSTATIONINFO__

// Latex documentation
#ifdef DOCSTRING_Loop
static const char *docstringLoopFileGnssStationInfo = R"(
\subsection{FileGnssStationInfo}
DEPRECATDED. Use LoopPlatformEquipment instead.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "classes/loop/loop.h"
#include "files/filePlatform.h"

/***** CLASS ***********************************/

/** @brief DEPRECATDED. Use LoopPlatformEquipment instead.
* @ingroup LoopGroup
* @see Loop */
class LoopFileGnssStationInfo : public Loop
{
  Platform                platform;
  PlatformEquipment::Type type;
  std::string             nameName, nameSerial, nameInfo, nameTimeStart, nameTimeEnd;
  std::string             nameIndex, nameCount;
  UInt                    index;

public:
  LoopFileGnssStationInfo(Config &config);

  UInt count() const override {return std::count_if(platform.equipments.begin(), platform.equipments.end(), [&](const auto &x)
                                                   {return (x->getType() == type);});}
  Bool iteration(VariableList &varList) override;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline LoopFileGnssStationInfo::LoopFileGnssStationInfo(Config &config)
{
  try
  {
    FileName    fileName;
    std::string choice;

    readConfig(config, "inputfileGnssStationInfo", fileName, Config::MUSTSET, "", "station/transmitter info file");
    if(readConfigChoice(config, "infoType", choice, Config::MUSTSET, "", "info to loop over"))
    {
      if(readConfigChoiceElement(config, "antenna",  choice, "loop over antennas"))  type = PlatformEquipment::GNSSANTENNA;
      if(readConfigChoiceElement(config, "receiver", choice, "loop over receivers")) type = PlatformEquipment::GNSSRECEIVER;
      endChoice(config);
    }
    readConfig(config, "variableLoopName",      nameName,      Config::OPTIONAL,  "loopName",      "variable with antenna/receiver name");
    readConfig(config, "variableLoopSerial",    nameSerial,    Config::OPTIONAL,  "loopSerial",    "variable with antenna/receiver serial");
    readConfig(config, "variableLoopInfo",      nameInfo,      Config::OPTIONAL,  "loopInfo",      "variable with radome (antenna) or version (receiver)");
    readConfig(config, "variableLoopTimeStart", nameTimeStart, Config::OPTIONAL,  "loopTimeStart", "variable with antenna/receiver start time");
    readConfig(config, "variableLoopTimeEnd",   nameTimeEnd,   Config::OPTIONAL,  "loopTimeEnd",   "variable with antenna/receiver end time");
    readConfig(config, "variableLoopIndex",     nameIndex,     Config::OPTIONAL,  "",              "variable with index of current iteration (starts with zero)");
    readConfig(config, "variableLoopCount",     nameCount,     Config::OPTIONAL,  "",              "variable with total number of iterations");
    if(isCreateSchema(config)) return;

    logWarningOnce<<"LoopFileGnssStationInfo is DEPRECATDED. Use LoopPlatformEquipment instead."<<Log::endl;

    readFilePlatform(fileName, platform);
    index = 0;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Bool LoopFileGnssStationInfo::iteration(VariableList &varList)
{
  if(index >= count())
    return FALSE;

  for(const auto &eq : platform.equipments)
    if(eq->getType() == type)
    {
      if(!nameIndex.empty())     addVariable(nameIndex,     index,               varList);
      if(!nameCount.empty())     addVariable(nameCount,     count(),             varList);
      if(!nameName.empty())      addVariable(nameName,      eq->name,            varList);
      if(!nameSerial.empty())    addVariable(nameSerial,    eq->serial,          varList);
      if(!nameTimeStart.empty()) addVariable(nameTimeStart, eq->timeStart.mjd(), varList);
      if(!nameTimeEnd.empty())   addVariable(nameTimeEnd,   eq->timeEnd.mjd(),   varList);
      if(!nameInfo.empty())
      {
        switch(eq->getType())
        {
          case PlatformEquipment::GNSSANTENNA:  addVariable(nameInfo, std::dynamic_pointer_cast<PlatformGnssAntenna>(eq)->radome,   varList); break;
          case PlatformEquipment::GNSSRECEIVER: addVariable(nameInfo, std::dynamic_pointer_cast<PlatformGnssReceiver>(eq)->version, varList); break;
          case PlatformEquipment::OTHER:        break;
          case PlatformEquipment::UNDEFINED:    break;
          default:                              break;
        }
      }
    }

  index++;
  return TRUE;
}

/***********************************************/

#endif
