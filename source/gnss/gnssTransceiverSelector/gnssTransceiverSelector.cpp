/***********************************************/
/**
* @file gnssTransceiverSelector.cpp
*
* @brief Selected receivers or transmitters.
*
* @author Torsten Mayer-Guerr
* @date 2021-01-23
*
*/
/***********************************************/

#define DOCSTRING_GnssTransceiverSelector

#include "base/import.h"
#include "config/configRegister.h"
#include "gnss/gnssTransceiverSelector/gnssTransceiverSelectorAll.h"
#include "gnss/gnssTransceiverSelector/gnssTransceiverSelectorExclude.h"
#include "gnss/gnssTransceiverSelector/gnssTransceiverSelectorFile.h"
#include "gnss/gnssTransceiverSelector/gnssTransceiverSelectorWildcard.h"
#include "gnss/gnssTransceiverSelector/gnssTransceiverSelector.h"

/***********************************************/

GROOPS_REGISTER_CLASS(GnssTransceiverSelector, "gnssTransceiverSelectorType",
                      GnssTransceiverSelectorAll,
                      GnssTransceiverSelectorWildcard,
                      GnssTransceiverSelectorFile,
                      GnssTransceiverSelectorExclude)

GROOPS_READCONFIG_UNBOUNDED_CLASS(GnssTransceiverSelector, "gnssTransceiverSelectorType")

/***********************************************/

GnssTransceiverSelector::GnssTransceiverSelector(Config &config, const std::string &name)
{
  try
  {
    std::string choice;
    while(readConfigChoice(config, name, choice, Config::OPTIONAL, "", "selected receivers or transmitters"))
    {
      if(readConfigChoiceElement(config, "all",      choice, "all available receivers/transmitters"))
        bases.push_back(std::unique_ptr<GnssTransceiverSelectorBase>(new GnssTransceiverSelectorAll(config)));
      if(readConfigChoiceElement(config, "wildcard", choice, "select by name and number"))
        bases.push_back(std::unique_ptr<GnssTransceiverSelectorBase>(new GnssTransceiverSelectorWildcard(config)));
      if(readConfigChoiceElement(config, "file",     choice, "select from file (with alternatives)"))
        bases.push_back(std::unique_ptr<GnssTransceiverSelectorBase>(new GnssTransceiverSelectorFile(config)));
      if(readConfigChoiceElement(config, "exclude",  choice, "exclude from selection"))
        bases.push_back(std::unique_ptr<GnssTransceiverSelectorBase>(new GnssTransceiverSelectorExclude(config)));
      endChoice(config);
      if(isCreateSchema(config))
        return;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::vector<Byte> GnssTransceiverSelector::select(const std::vector<GnssTransceiverPtr> &transceivers) const
{
  try
  {
    std::vector<Byte> selected(transceivers.size(), (bases.size() && bases.front()->exclude()));
    for(auto &base : bases)
      base->select(transceivers, selected);
    for(UInt i=0; i<transceivers.size(); i++)
      selected.at(i) = selected.at(i) && transceivers.at(i)->useable();
    return selected;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::vector<Byte> GnssTransceiverSelector::select(const std::vector<GnssReceiverPtr> &receivers) const
{
  try
  {
    std::vector<GnssTransceiverPtr> transceivers;
    transceivers.insert(transceivers.begin(), receivers.begin(), receivers.end());
    return select(transceivers);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::vector<Byte> GnssTransceiverSelector::select(const std::vector<GnssTransmitterPtr> &transmitters) const
{
  try
  {
    std::vector<GnssTransceiverPtr> transceivers;
    transceivers.insert(transceivers.begin(), transmitters.begin(), transmitters.end());
    return select(transceivers);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
