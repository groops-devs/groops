/***********************************************/
/**
* @file gnssParametrizationReceiver.h
*
* @brief GNSS receivers.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2020-06-02
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPARAMETRIZATIONRECEIVER__
#define __GROOPS_GNSSPARAMETRIZATIONRECEIVER__

// Latex documentation
#ifdef DOCSTRING_GnssParametrizationReceiver
static const char *docstringGnssParametrizationReceiver = R"(
\section{GnssParametrizationReceiver}\label{gnssParametrizationReceiverType}

Definition and parametrization of GNSS receivers.

Input files used in \config{stationNetwork} and \config{lowEarthOrbiter}:
\begin{itemize}
  \item \configFile{inputfileStationInfo}{gnssStationInfo}:
        Created via \program{GnssStationLog2StationInfo} or \program{GnssStationInfoCreate}.
  \item \configFile{inputfileAntennaDefinition}{gnssAntennaDefinition}:
        Created via \program{GnssAntex2AntennaDefinition} or \program{GnssAntennaDefinitionCreate}.
  \item \configFile{inputfileReceiverDefinition}{gnssReceiverDefinition}:
        Created via \program{GnssReceiverDefinitionCreate} in case you want to define which signal
        types a receiver model can observe.
  \item \configFile{inputfileAccuracyDefinition}{gnssAntennaDefinition}:
        Created via \program{GnssAntennaDefinitionCreate}.
  \item \configFile{inputfileSignalBias}{gnssSignalBias}:
        Estimated via \program{GnssProcessing}.
  \item \configFile{inputfileObservation}{instrument}:
        Converted from RINEX observation files via \program{RinexObservation2GnssReceiver}.
\end{itemize}

See also \program{GnssProcessing} and \program{GnssSimulateReceiver}.
)";
#endif

/***********************************************/

#include "config/config.h"
#include "gnss/gnss.h"

/**
* @defgroup gnssParametrizationReceiverGroup GnssParametrizationReceiver
* @brief GNSS receivers.
* @ingroup gnssGroup
* Parametrization of GNSS receivers.
* The interface is given by @ref GnssParametrizationReceiver.
* An Instance can be created by @ref readConfig. */
/// @{

/***** TYPES ***********************************/

class GnssParametrizationReceiver;
typedef std::shared_ptr<GnssParametrizationReceiver> GnssParametrizationReceiverPtr;

/***** CLASS ***********************************/

class GnssParametrizationReceiver
{
public:
  virtual ~GnssParametrizationReceiver() {}

  /** @brief List of receivers. */
  virtual std::vector<Gnss::ReceiverPtr> receivers() = 0;

  /** @brief List of parametrizations. */
  virtual std::vector<Gnss::ParametrizationPtr> parametrizations() = 0;

  /** @brief creates an derived instance of this class. */
  static GnssParametrizationReceiverPtr create(Config &config, const std::string &name);
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class GnssParametrizationReceiver.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a var is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] var Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates GnssParametrizationReceiver */
template<> Bool readConfig(Config &config, const std::string &name, GnssParametrizationReceiverPtr &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***********************************************/

#endif /* __GROOPS___ */
