/***********************************************/
/**
* @file platformSelector.cpp
*
* @brief Selected platforms.
*
* @author Torsten Mayer-Guerr
* @date 2021-01-23
*
*/
/***********************************************/

#define DOCSTRING_PlatformSelector

#include "base/import.h"
#include "config/configRegister.h"
#include "classes/platformSelector/platformSelectorAll.h"
#include "classes/platformSelector/platformSelectorExclude.h"
#include "classes/platformSelector/platformSelectorFile.h"
#include "classes/platformSelector/platformSelectorWildcard.h"
#include "classes/platformSelector/platformSelector.h"

/***********************************************/

GROOPS_REGISTER_CLASS(PlatformSelector, "platformSelectorType",
                      PlatformSelectorAll,
                      PlatformSelectorWildcard,
                      PlatformSelectorFile,
                      PlatformSelectorExclude)

GROOPS_READCONFIG_UNBOUNDED_CLASS(PlatformSelector, "platformSelectorType")

/***********************************************/

PlatformSelector::PlatformSelector(Config &config, const std::string &name)
{
  try
  {
    std::string choice;
    while(readConfigChoice(config, name, choice, Config::OPTIONAL, "", "selected platforms (stations, satellites, ...)"))
    {
      if(readConfigChoiceElement(config, "all",      choice, "all available platforms"))
        bases.push_back(std::unique_ptr<PlatformSelectorBase>(new PlatformSelectorAll(config)));
      if(readConfigChoiceElement(config, "wildcard", choice, "select by name and number"))
        bases.push_back(std::unique_ptr<PlatformSelectorBase>(new PlatformSelectorWildcard(config)));
      if(readConfigChoiceElement(config, "file",     choice, "select from file (with alternatives)"))
        bases.push_back(std::unique_ptr<PlatformSelectorBase>(new PlatformSelectorFile(config)));
      if(readConfigChoiceElement(config, "exclude",  choice, "exclude from selection"))
        bases.push_back(std::unique_ptr<PlatformSelectorBase>(new PlatformSelectorExclude(config)));
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

std::vector<Byte> PlatformSelector::select(const Time &timeStart, const Time &timeEnd, const std::vector<const Platform*> &platforms) const

{
  try
  {
    std::vector<Byte> selected(platforms.size(), (bases.size() && bases.front()->exclude()));
    for(auto &base : bases)
      base->select(timeStart, timeEnd, platforms, selected);
    for(UInt i=0; i<platforms.size(); i++)
      selected.at(i) = selected.at(i) && platforms.at(i);
    return selected;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
