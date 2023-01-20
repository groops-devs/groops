/***********************************************/
/**
* @file platformSelectorWildcard.h
*
* @brief Select platforms with wildcards.
* @see PlatformSelector
*
* @author Torsten Mayer-Guerr
* @date 2021-01-23
*
*/
/***********************************************/

#ifndef __GROOPS_PLATFORMSELECTORWILDCARD__
#define __GROOPS_PLATFORMSELECTORWILDCARD__

// Latex documentation
#ifdef DOCSTRING_PlatformSelector
static const char *docstringPlatformSelectorWildcard = R"(
\subsection{Wildcard}\label{platformSelectorType:wildcard}
Select all receivers/transmitters which match the
\config{name}, \config{markerName}, and \config{markerNumber}.
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
class PlatformSelectorWildcard : public PlatformSelectorBase
{
  std::regex patternName, patternMarkerName, pattermMarkerNumber;

public:
  PlatformSelectorWildcard(Config &config);
  void select(const Time &timeStart, const Time &timeEnd, const std::vector<const Platform*> &platforms, std::vector<Byte> &selected) const override;
};

/***********************************************/

inline PlatformSelectorWildcard::PlatformSelectorWildcard(Config &config)
{
  try
  {
    std::string name, markerName, markerNumber;

    readConfig(config, "name",         name,         Config::OPTIONAL, "*", "wildcards: * and ?");
    readConfig(config, "markerName",   markerName,   Config::OPTIONAL, "*", "wildcards: * and ?, from platform");
    readConfig(config, "markerNumber", markerNumber, Config::OPTIONAL, "*", "wildcards: * and ?, from platform");
    if(isCreateSchema(config)) return;

    patternName         = String::wildcard2regex(name);
    patternMarkerName   = String::wildcard2regex(markerName);
    pattermMarkerNumber = String::wildcard2regex(markerNumber);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void PlatformSelectorWildcard::select(const Time &/*timeStart*/, const Time &/*timeEnd*/, const std::vector<const Platform*> &platforms, std::vector<Byte> &selected) const
{
  try
  {
    for(UInt i=0; i<platforms.size(); i++)
      if(platforms.at(i) &&
         std::regex_match(platforms.at(i)->name,         patternName) &&
         std::regex_match(platforms.at(i)->markerName,   patternMarkerName) &&
         std::regex_match(platforms.at(i)->markerNumber, pattermMarkerNumber))
        selected.at(i) = TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
