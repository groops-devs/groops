/***********************************************/
/**
* @file gnssTransmitterGenerator.h
*
* @brief Provides a list of transmitters.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2021-02-25
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSTRANSMITTERGENERATOR__
#define __GROOPS_GNSSTRANSMITTERGENERATOR__

// Latex documentation
#ifdef DOCSTRING_GnssTransmitterGenerator
static const char *docstringGnssTransmitterGenerator = R"(
\section{GnssTransmitterGenerator}\label{gnssTransmitterGeneratorType}
Definition and basic information of GNSS transmitters.

See also \program{GnssProcessing} and \program{GnssSimulateReceiver}.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "gnss/gnssTransmitter.h"

/**
* @defgroup gnssTransmitterGeneratorGroup GnssTransmitterGenerator
* @brief Provides a list of transmitters.
* @ingroup classesGroup
* The interface is given by @ref GnssTransmitterGenerator.
* An Instance can be created by @ref readConfig. */
/// @{

/***** TYPES ***********************************/

class GnssTransmitterGenerator;
class GnssTransmitterGeneratorBase;
typedef std::shared_ptr<GnssTransmitterGenerator> GnssTransmitterGeneratorPtr;

/***** CLASS ***********************************/

/** @brief Provides a list of transmitters.
* An Instance of this class can be created by @ref readConfig. */
class GnssTransmitterGenerator
{
  std::vector<GnssTransmitterGeneratorBase*> base;

public:
  /** @brief Constructor from config. */
  GnssTransmitterGenerator(Config &config, const std::string &name);

  /// Destructor.
 ~GnssTransmitterGenerator();

  /** @brief Initialize and returns a vector of transmitters. */
  std::vector<GnssTransmitterPtr> transmitters(const std::vector<Time> &times);

  /** @brief creates an derived instance of this class. */
  static GnssTransmitterGeneratorPtr create(Config &config, const std::string &name) {return GnssTransmitterGeneratorPtr(new GnssTransmitterGenerator(config, name));}
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class GnssTransmitterGenerator.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a var is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] var Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates GnssTransmitterGenerator */
template<> Bool readConfig(Config &config, const std::string &name, GnssTransmitterGeneratorPtr &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***** CLASS ***********************************/

// Internal class
class GnssTransmitterGeneratorBase
{
public:
  virtual ~GnssTransmitterGeneratorBase() {}
  virtual void init(const std::vector<Time> &times, std::vector<GnssTransmitterPtr> &transmitters) = 0;
};

/***********************************************/

#endif /* __GROOPS___ */
