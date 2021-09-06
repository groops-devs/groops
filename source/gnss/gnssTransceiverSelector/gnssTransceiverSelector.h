/***********************************************/
/**
* @file gnssTransceiverSelector.h
*
* @brief Selected receivers or transmitters.
*
* @author Torsten Mayer-Guerr
* @date 2021-01-23
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSTRANSCEIVERSELECTOR__
#define __GROOPS_GNSSTRANSCEIVERSELECTOR__

// Latex documentation
#ifdef DOCSTRING_GnssTransceiverSelector
static const char *docstringGnssTransceiverSelector = R"(
\section{GnssTransceiverSelector}\label{gnssTransceiverSelectorType}
Select a list of GNSS transmitters or receivers.

See also \program{GnssProcessing}.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "gnss/gnssReceiver.h"
#include "gnss/gnssTransmitter.h"

/**
* @defgroup gnssTransceiverSelectorGroup GnssTransceiverSelector
* @brief Selected receivers or transmitters.
* @ingroup classesGroup
* The interface is given by @ref GnssTransceiverSelector.
* An Instance can be created by @ref readConfig. */
/// @{

/***** TYPES ***********************************/

class GnssTransceiverSelector;
class GnssTransceiverSelectorBase;
typedef std::shared_ptr<GnssTransceiverSelector> GnssTransceiverSelectorPtr;

/***** CLASS ***********************************/

/** @brief Selected receivers or transmitters.
* An instance of this class can be created with @ref readConfig. */
class GnssTransceiverSelector
{
  std::vector<std::unique_ptr<GnssTransceiverSelectorBase>> bases;

public:
  /// Constructor.
  GnssTransceiverSelector(Config &config, const std::string &name);

  /** @brief returns a boolean vector which receivers/transmitters are selected. */
  std::vector<Byte> select(const std::vector<GnssTransceiverPtr> &transceivers) const;
  std::vector<Byte> select(const std::vector<GnssReceiverPtr>    &receivers)    const;
  std::vector<Byte> select(const std::vector<GnssTransmitterPtr> &transmitters) const;

  /** @brief creates an derived instance of this class. */
  static GnssTransceiverSelectorPtr create(Config &config, const std::string &name) {return GnssTransceiverSelectorPtr(new GnssTransceiverSelector(config, name));}
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class GnssTransceiverSelector.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and an class without points is created.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] var Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates GnssTransceiverSelector */
template<> Bool readConfig(Config &config, const std::string &name, GnssTransceiverSelectorPtr &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***** CLASS ***********************************/

// Internal class
class GnssTransceiverSelectorBase
{
public:
  virtual ~GnssTransceiverSelectorBase() {}
  virtual void select(const std::vector<GnssTransceiverPtr> &transceivers, std::vector<Byte> &selected) const = 0;
  virtual Bool exclude() const {return FALSE;}
};

/***********************************************/

#endif
