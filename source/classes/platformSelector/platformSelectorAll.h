/***********************************************/
/**
* @file platformSelectorAll.h
*
* @brief All available platforms.
* @see PlatformSelector
*
* @author Torsten Mayer-Guerr
* @date 2021-01-23
*
*/
/***********************************************/

#ifndef __GROOPS_PLATFORMSELECTORALL__
#define __GROOPS_PLATFORMSELECTORALL__

// Latex documentation
#ifdef DOCSTRING_PlatformSelector
static const char *docstringPlatformSelectorAll = R"(
\subsection{All}\label{platformSelectorType:all}
Select all platforms.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "classes/platformSelector/platformSelector.h"

/***** CLASS ***********************************/

/** @brief All available platforms.
* @ingroup platformSelectorGroup
* @see PlatformSelector */
class PlatformSelectorAll : public PlatformSelectorBase
{
public:
  PlatformSelectorAll(Config &/*config*/) {}

  void select(const Time &timeStart, const Time &timeEnd, const std::vector<const Platform*> &platforms, std::vector<Byte> &selected) const override;
};

/***********************************************/

inline void PlatformSelectorAll::select(const Time &/*timeStart*/, const Time &/*timeEnd*/, const std::vector<const Platform*> &/*platforms*/, std::vector<Byte> &selected) const
{
  std::fill(selected.begin(), selected.end(), TRUE);
}

/***********************************************/

#endif
