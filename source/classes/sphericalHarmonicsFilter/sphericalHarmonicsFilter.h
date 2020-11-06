/***********************************************/
/**
* @file sphericalHarmonicsFilter.h
*
* @brief Filtering spherical harmonics.
*
* @author Torsten Mayer-Guerr
* @date 2008-06-08
*
*/
/***********************************************/

#ifndef __GROOPS_SPHERICALHARMONICSFILTER__
#define __GROOPS_SPHERICALHARMONICSFILTER__

// Latex documentation
#ifdef DOCSTRING_SphericalHarmonicsFilter
static const char *docstringSphericalHarmonicsFilter = R"(
\section{SphericalHarmonicsFilter}\label{sphericalHarmonicsFilterType}
Filtering of a spherical harmonics expansion.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "base/sphericalHarmonics.h"
#include "config/config.h"

/**
* @defgroup sphericalHarmonicsFilterGroup SphericalHarmonicsFilter
* @brief Filtering spherical harmonics.
* @ingroup classesGroup
* The interface is given by @ref SphericalHarmonicsFilter.
* An Instance can be created by @ref readConfig. */
/// @{

/***** TYPES ***********************************/

class SphericalHarmonicsFilter;
class SphericalHarmonicsFilterBase;
typedef std::shared_ptr<SphericalHarmonicsFilter> SphericalHarmonicsFilterPtr;

/***** CLASS ***********************************/

/** @brief Filtering spherical harmonics.
* An Instance of this class can be created by @ref readConfig. */
class SphericalHarmonicsFilter
{
  std::vector<SphericalHarmonicsFilterBase*> filters;

public:
  /// Constructor.
  SphericalHarmonicsFilter(Config &config, const std::string &name);

  /// Destructor.
  virtual ~SphericalHarmonicsFilter();

  /** @brief returns a filtered version of given harmonics. */
  virtual SphericalHarmonics filter(const SphericalHarmonics &harm) const;

  /** @brief creates an derived instance of this class. */
  static SphericalHarmonicsFilterPtr create(Config &config, const std::string &name) {return SphericalHarmonicsFilterPtr(new SphericalHarmonicsFilter(config, name));}
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class SphericalHarmonicsFilter.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and a filter with no effect is created.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] filter Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates SphericalHarmonicsFilter */
template<> Bool readConfig(Config &config, const std::string &name, SphericalHarmonicsFilterPtr &filter, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***** CLASS ***********************************/

// Internal class
class SphericalHarmonicsFilterBase
{
public:
virtual ~SphericalHarmonicsFilterBase() {}
virtual SphericalHarmonics filter(const SphericalHarmonics &harm) const = 0;
};

/***********************************************/

#endif /* __GROOPS_FILTER__ */
