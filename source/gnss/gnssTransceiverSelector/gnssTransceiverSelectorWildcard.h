/***********************************************/
/**
* @file gnssTransceiverSelectorWildcard.h
*
* @brief Select transceivers with wildcards.
* @see GnssTransceiverSelector
*
* @author Torsten Mayer-Guerr
* @date 2021-01-23
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSTRANSCEIVERSELECTORWILDCARD__
#define __GROOPS_GNSSTRANSCEIVERSELECTORWILDCARD__

// Latex documentation
#ifdef DOCSTRING_GnssTransceiverSelector
static const char *docstringGnssTransceiverSelectorWildcard = R"(
\subsection{Wildcard}\label{gnssTransceiverSelectorType:wildcard}
Select all receivers/transmitters which match the
\config{name}, \config{markerName}, and \config{markerNumber}.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "base/string.h"
#include "config/config.h"
#include "gnss/gnssTransceiverSelector/gnssTransceiverSelector.h"

/***** CLASS ***********************************/

/** @brief Select transceivers with wildcards.
* @ingroup gnssTransceiverSelectorGroup
* @see GnssTransceiverSelector */
class GnssTransceiverSelectorWildcard : public GnssTransceiverSelectorBase
{
  std::regex patternName, patternMarkerName, pattermMarkerNumber;

public:
  GnssTransceiverSelectorWildcard(Config &config);
  void select(const std::vector<GnssTransceiverPtr> &transceivers, std::vector<Byte> &selected) const override;
};

/***********************************************/

inline GnssTransceiverSelectorWildcard::GnssTransceiverSelectorWildcard(Config &config)
{
  try
  {
    std::string name, markerName, markerNumber;

    readConfig(config, "name",         name,         Config::OPTIONAL, "*", "wildcards: * and ?");
    readConfig(config, "markerName",   markerName,   Config::OPTIONAL, "*", "wildcards: * and ?, from stationInfo");
    readConfig(config, "markerNumber", markerNumber, Config::OPTIONAL, "*", "wildcards: * and ?, from stationInfo");
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

inline void GnssTransceiverSelectorWildcard::select(const std::vector<GnssTransceiverPtr> &transceivers, std::vector<Byte> &selected) const
{
  try
  {
    for(UInt i=0; i<transceivers.size(); i++)
      if(std::regex_match(transceivers.at(i)->name(),            patternName) &&
         std::regex_match(transceivers.at(i)->info.markerName,   patternMarkerName) &&
         std::regex_match(transceivers.at(i)->info.markerNumber, pattermMarkerNumber))
        selected.at(i) = TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
