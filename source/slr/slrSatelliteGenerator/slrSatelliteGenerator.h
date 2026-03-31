/***********************************************/
/**
* @file slrSatelliteGenerator.h
*
* @brief Provides a list of satellites.
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#ifndef __GROOPS_SLRSATELLITEGENERATOR__
#define __GROOPS_SLRSATELLITEGENERATOR__

// Latex documentation
#ifdef DOCSTRING_SlrSatelliteGenerator
static const char *docstringSlrSatelliteGenerator = R"(
\section{SlrSatelliteGenerator}\label{slrSatelliteGeneratorType}
Definition and basic information of SLR satellites.

See also \program{SlrProcessing}.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "slr/slrSatellite.h"

/**
* @defgroup slrSatelliteGeneratorGroup SlrSatelliteGenerator
* @brief Provides a list of satellites.
* @ingroup classesGroup
* The interface is given by @ref SlrSatelliteGenerator. */
/// @{

/***** TYPES ***********************************/

class SlrSatelliteGenerator;
class SlrSatelliteGeneratorBase;
typedef std::shared_ptr<SlrSatelliteGenerator> SlrSatelliteGeneratorPtr;

/***** CLASS ***********************************/

/** @brief Provides a list of satellites.
* An Instance of this class can be created by @ref readConfig. */
class SlrSatelliteGenerator
{
  std::vector<SlrSatelliteGeneratorBase*> base;

public:
  /** @brief Constructor from config. */
  SlrSatelliteGenerator(Config &config, const std::string &name);

  /// Destructor.
 ~SlrSatelliteGenerator();

  /** @brief Iniatialize and returns a vector of satellites. */
  std::vector<SlrSatellitePtr> satellites(const std::vector<Time> &times);

  /** @brief creates an derived instance of this class. */
  static SlrSatelliteGeneratorPtr create(Config &config, const std::string &name) {return SlrSatelliteGeneratorPtr(new SlrSatelliteGenerator(config, name));}
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class SlrSatelliteGenerator.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a var is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] var Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates SlrSatelliteGenerator */
template<> Bool readConfig(Config &config, const std::string &name, SlrSatelliteGeneratorPtr &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***** CLASS ***********************************/

// Internal class
class SlrSatelliteGeneratorBase
{
public:
  virtual ~SlrSatelliteGeneratorBase() {}
  virtual void init(const std::vector<Time> &times, std::vector<SlrSatellitePtr> &satellites) = 0;
};

/***********************************************/

#endif /* __GROOPS___ */
