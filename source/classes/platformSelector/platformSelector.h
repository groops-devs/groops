/***********************************************/
/**
* @file platformSelector.h
*
* @brief Selected platforms.
*
* @author Torsten Mayer-Guerr
* @date 2021-01-23
*
*/
/***********************************************/

#ifndef __GROOPS_PLATFORMSELECTOR__
#define __GROOPS_PLATFORMSELECTOR__

// Latex documentation
#ifdef DOCSTRING_PlatformSelector
static const char *docstringPlatformSelector = R"(
\section{PlatformSelector}\label{platformSelectorType}
Select a list of platforms (stations, satellites, ...).
In a first step all platforms are selected if first selector \config{exclude}s platforms
otherwise all platforms excluded. When every selector from top to bottom selects or deselects
(with \config{exclude}) the matching platforms.

See also \program{GnssProcessing} or \program{SlrProcessing}.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "files/filePlatform.h"

/**
* @defgroup platformSelectorGroup PlatformSelector
* @brief Selected platforms.
* @ingroup classesGroup
* The interface is given by @ref PlatformSelector. */
/// @{

/***** TYPES ***********************************/

class PlatformSelector;
class PlatformSelectorBase;
typedef std::shared_ptr<PlatformSelector> PlatformSelectorPtr;

/***** CLASS ***********************************/

/** @brief Selected platforms.
* An instance of this class can be created with @ref readConfig. */
class PlatformSelector
{
  std::vector<std::unique_ptr<PlatformSelectorBase>> bases;

public:
  /// Constructor.
  PlatformSelector(Config &config, const std::string &name);

  /** @brief returns a boolean vector which platforms are selected. */
  std::vector<Byte> select(const Time &timeStart, const Time &timeEnd, const std::vector<const Platform*> &platforms) const;

  /** @brief creates an derived instance of this class. */
  static PlatformSelectorPtr create(Config &config, const std::string &name) {return PlatformSelectorPtr(new PlatformSelector(config, name));}
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class PlatformSelector.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and an class without points is created.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] var Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates PlatformSelector */
template<> Bool readConfig(Config &config, const std::string &name, PlatformSelectorPtr &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***** CLASS ***********************************/

// Internal class
class PlatformSelectorBase
{
public:
  Bool exclude;

  virtual ~PlatformSelectorBase() {}
  virtual void select(const Time &timeStart, const Time &timeEnd, const std::vector<const Platform*> &platforms, std::vector<Byte> &selected) const = 0;
};

/***********************************************/

#endif
