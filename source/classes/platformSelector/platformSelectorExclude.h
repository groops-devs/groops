/***********************************************/
/**
* @file platformSelectorExclude.h
*
* @brief Exclude platforms.
* @see PlatformSelector
*
* @author Torsten Mayer-Guerr
* @date 2021-01-23
*
*/
/***********************************************/

#ifndef __GROOPS_PLATFORMSELECTOREXCLUDE__
#define __GROOPS_PLATFORMSELECTOREXCLUDE__

// Latex documentation
#ifdef DOCSTRING_PlatformSelector
static const char *docstringPlatformSelectorExclude = R"(
\subsection{Exclude}\label{platformSelectorType:exclude}
Deselects all selected receivers/transmitters of
\configClass{selector}{platformSelectorType}.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "classes/platformSelector/platformSelector.h"

/***** CLASS ***********************************/

/** @brief Exclude platforms.
* @ingroup platformSelectorGroup
* @see PlatformSelector */
class PlatformSelectorExclude : public PlatformSelectorBase
{
  PlatformSelectorPtr selector;

public:
  PlatformSelectorExclude(Config &config);
  void select(const Time &timeStart, const Time &timeEnd, const std::vector<const Platform*> &platforms, std::vector<Byte> &selected) const override;
};

/***********************************************/

inline PlatformSelectorExclude::PlatformSelectorExclude(Config &config)
{
  try
  {
    readConfig(config, "selector", selector, Config::MUSTSET, "", "");
    exclude = TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void PlatformSelectorExclude::select(const Time &timeStart, const Time &timeEnd, const std::vector<const Platform*> &platforms, std::vector<Byte> &selected) const
{
  try
  {
    std::vector<Byte> excl = selector->select(timeStart, timeEnd, platforms);
    for(UInt i=0; i<excl.size(); i++)
      if(excl.at(i))
        selected.at(i) = FALSE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
