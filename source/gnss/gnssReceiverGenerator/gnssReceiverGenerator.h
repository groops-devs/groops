/***********************************************/
/**
* @file gnssReceiverGenerator.h
*
* @brief Provides a list of receivers.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2021-02-25
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSRECEIVERGENERATOR__
#define __GROOPS_GNSSRECEIVERGENERATOR__

// Latex documentation
#ifdef DOCSTRING_GnssReceiverGenerator
static const char *docstringGnssReceiverGenerator = R"(
\section{GnssReceiverGenerator}\label{gnssReceiverGeneratorType}
Definition and basic information of GNSS receivers.

Most of the input files are provided in GROOPS file formats at
\url{https://ftp.tugraz.at/outgoing/ITSG/groops} (marked with \textbf{*} below).
These files are regularly updated.
\begin{itemize}
  \item \configFile{inputfileStationInfo}{gnssStationInfo}\textbf{*}:
        Antenna and receiver information, antenna reference point offsets, antenna orientations.
        Created via \program{GnssStationLog2StationInfo} or \program{GnssStationInfoCreate}.
  \item \configFile{inputfileAntennaDefinition}{gnssAntennaDefinition}\textbf{*}:
        Antenna center offsets and variations.
        Created via \program{GnssAntex2AntennaDefinition} or \program{GnssAntennaDefinitionCreate}.
  \item \configFile{inputfileReceiverDefinition}{gnssReceiverDefinition}:
        Observed signal types (optional).
        Created via \program{GnssReceiverDefinitionCreate} in case you want to define which signal
        types a receiver model can observe.
  \item \configFile{inputfileAccuracyDefinition}{gnssAntennaDefinition}\textbf{*}:
        Elevation and azimuth dependent accuracy.
        Created via \program{GnssAntennaDefinitionCreate}.
  \item \configFile{inputfileObservation}{instrument}:
        Converted from RINEX observation files via \program{RinexObservation2GnssReceiver}.
\end{itemize}

It is possible to limit the observation types to be used in the processing by a list of \configClass{useType}{gnssType}
and any observation types not defined within the list are ignored and discarded.
Similarly observations defined in the list of \configClass{ignoreType}{gnssType} are ignored and discarded.
The codes used follow the \href{https://files.igs.org/pub/data/format/rinex305.pdf}{RINEX 3 definition}.

Each receiver goes through a \config{preprocessing} step individually, where observation outliers are removed or downweighted,
continuous tracks of phase observations are defined for ambiguity parametrization, cycle slips are detected, and receivers are
disabled if they do not fulfill certain requirements. The preprocessing step consists of an initial PPP estimation done by
\reference{robust least squares adjustment}{fundamentals.robustLeastSquares} and checks whether the position error
of the solutions exceeds \config{codeMaxPositionDiff}. If the error exceeds the threshold the receiver will be discarded.
The preprocessing also sets initial clock error values and removes tracks that stay below a certain elevation mask (\config{elevationTrackMinimum}).

See also \program{GnssProcessing} and \program{GnssSimulateReceiver}.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/noiseGenerator/noiseGenerator.h"
#include "gnss/gnssReceiver.h"
#include "gnss/gnssTransmitter.h"

/**
* @defgroup gnssReceiverGeneratorGroup GnssReceiverGenerator
* @brief Provides a list of receivers.
* @ingroup classesGroup
* The interface is given by @ref GnssReceiverGenerator.
* An Instance can be created by @ref readConfig. */
/// @{

/***** TYPES ***********************************/

class Gnss;
class GnssReceiverGenerator;
class GnssReceiverGeneratorBase;
typedef std::shared_ptr<GnssReceiverGenerator> GnssReceiverGeneratorPtr;

/***** CLASS ***********************************/

/** @brief Provides a list of receivers.
* An Instance of this class can be created by @ref readConfig. */
class GnssReceiverGenerator
{
  std::vector<GnssReceiverGeneratorBase*> base;

public:
  /** @brief Constructor from config. */
  GnssReceiverGenerator(Config &config, const std::string &name);

  /// Destructor.
 ~GnssReceiverGenerator();

  /** @brief Iniatialize and returns a vector of receivers. */
  std::vector<GnssReceiverPtr> receivers(const std::vector<Time> &times, const Time &timeMargin,
                                         const std::vector<GnssTransmitterPtr> &transmitters, EarthRotationPtr earthRotation,
                                         Parallel::CommunicatorPtr comm);

  /** @brief preprocess the observations of receivers. */
  void preprocessing(Gnss *gnss, Parallel::CommunicatorPtr comm);

  /** @brief simulate the observations of receivers. */
  void simulation(const std::vector<GnssType> &types, NoiseGeneratorPtr noiseClock, NoiseGeneratorPtr noiseObs,
                  Gnss *gnss, Parallel::CommunicatorPtr comm);

  /** @brief creates an derived instance of this class. */
  static GnssReceiverGeneratorPtr create(Config &config, const std::string &name) {return GnssReceiverGeneratorPtr(new GnssReceiverGenerator(config, name));}
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class GnssReceiverGenerator.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a var is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] var Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates GnssReceiverGenerator */
template<> Bool readConfig(Config &config, const std::string &name, GnssReceiverGeneratorPtr &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***** CLASS ***********************************/

// Internal class
class GnssReceiverGeneratorBase
{
public:
  virtual ~GnssReceiverGeneratorBase() {}

  virtual void init(const std::vector<Time> &times, const Time &timeMargin, const std::vector<GnssTransmitterPtr> &transmitters,
                    EarthRotationPtr earthRotation, Parallel::CommunicatorPtr comm, std::vector<GnssReceiverPtr> &receivers) = 0;

  virtual void preprocessing(Gnss *gnss, Parallel::CommunicatorPtr comm) = 0;

  virtual void simulation(const std::vector<GnssType> &types, NoiseGeneratorPtr noiseClock, NoiseGeneratorPtr noiseObs,
                          Gnss *gnss, Parallel::CommunicatorPtr comm) = 0;

  static void printPreprocessingInfos(const std::string &header, const std::vector<GnssReceiverPtr> &receivers,
                                      Bool disabledOnly, Parallel::CommunicatorPtr comm);
};

/***********************************************/

#endif /* __GROOPS___ */
