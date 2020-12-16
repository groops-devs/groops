/***********************************************/
/**
* @file gnssParametrizationTransmitter.h
*
* @brief transmitter satellite constellations.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2018-08-07
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPARAMETRIZATIONTRANSMITTER__
#define __GROOPS_GNSSPARAMETRIZATIONTRANSMITTER__

// Latex documentation
#ifdef DOCSTRING_GnssParametrizationTransmitter
static const char *docstringGnssParametrizationTransmitter = R"(
\section{GnssParametrizationTransmitter}\label{gnssParametrizationTransmitterType}

Definition and parametrization of GNSS satelllite constellations. For precise point positioning (PPP) the
constellations are held fixed and no additional parameters are set up. In case of
\reference{GNSS satellite orbit determination}{cookbook.gnssNetwork} orbits, clock errors, and signal biases
are usually estimated.

A list of satellite PRNs (i.e for GPS: G01, G02, G03, ...) must be provided via
\configFile{inputfileTransmitterList}{stringList}. Satellite system codes follow the
\href{https://files.igs.org/pub/data/format/rinex304.pdf}{RINEX 3 definition}, see \reference{GnssType}{gnssType}.
All input and output files except \configFile{inputfileAntennaDefinition}{gnssAntennaDefinition},
\configFile{inputfileReceiverDefinition}{gnssReceiverDefinition}, and
\configFile{outputfileUsedTransmitterList}{stringList} are read/written for each satellite.
The file name is interpreted as a template with the variable \verb|{prn}| being replaced by the satellite PRN.

In case of PPP precise orbits must be provided via \configFile{inputfileVariational}{variationalEquation}.
Otherwise, dynamic orbits can be estimated by setting \config{estimateOrbit:dynamic}, with
\configFile{inputfileVariational}{variationalEquation} coming from
\reference{orbit integration}{cookbook.gnssNetwork:orbitIntegration}.
The \configClass{parametrizationAcceleration}{parametrizationAccelerationType} must include at least those
parameters that were estimated in \program{PreprocessingVariationalEquationOrbitFit}.
Additional \config{stochasticPulse} parameters can be set up to reduce orbit mismodeling.

Precise satellite clock errors have to be provided via \configFile{inputfileClock}{instrument} in case of PPP.
Otherwise, they can be estimated via \config{estimateClockError}, either \config{epochWise} or based
on a \config{noiseModel}. In this case it is enough to provide approximate clock errors (e.g. from broadcast ephemeris)
via \configFile{inputfileClock}{instrument}. The clock datum can be defined via a zero-mean constraint over all satellites
of a constellation using \config{sigmaZeroMeanConstraint}. This constraint is applied per epoch.

Signal biases for code and phase observations have to be provided via \configFile{inputfileSignalBias}{gnssSignalBias}
in case of PPP. Otherwise they can be estimated by setting \config{estimateCodeBias}=\verb|yes| and
\config{estimatePhaseBias}=\verb|yes|, in which case no a priori signal biases are needed.

Antenna center offsets and variations can be estimated by defining \configClass{antennaCenter}{gnssParametrizationAntennaType}.
It is possible to estimate patterns for each individual satellite or grouped by satellite block.
Since it is not possible to estimate, for example, antenna center offets based on a single day, they are usually hard-constrained
to their a priori value. This can be done by adding \configClass{constraints}{gnssParametrizationConstraintsType} in \program{GnssProcessing}.

Input files used in all systems:
\begin{itemize}
  \item \configFile{inputfileTransmitterInfo}{gnssStationInfo}:
        Created via \program{GnssAntex2AntennaDefinition} or \program{GnssStationInfoCreate}.
  \item \configFile{inputfileAntennaDefinition}{gnssAntennaDefinition}:
        Created via \program{GnssAntex2AntennaDefinition} or \program{GnssAntennaDefinitionCreate}.
  \item \configFile{inputfileReceiverDefinition}{gnssReceiverDefinition}:
        Created via \program{GnssReceiverDefinitionCreate} in case you want to define which signal
        types a satellite transmits.
  \item \configFile{inputfileSignalBias}{gnssSignalBias}: Converted via (TODO: SinexBias2GnssSignalBias) or estimated via \program{GnssProcessing}.
  \item \configFile{inputfileVariational}{variationalEquation}:
        Integrated via \program{PreprocessingVariationalEquation} and \program{PreprocessingVariationalEquationOrbitFit},
        created via (TODO: GnssVariationalCreate), or output of \program{GnssProcessing}.
  \item \configFile{inputfileClock}{instrument}:
        Converted via \program{GnssClockRinex2InstrumentClock} or \program{GnssRinexNavigation2OrbitClock} or
        output of \program{GnssProcessing}.
\end{itemize}

See also \program{GnssProcessing} and \program{GnssSimulateReceiver}.
)";
#endif

/***********************************************/

#include "base/polynomial.h"
#include "files/fileGnssStationInfo.h"
#include "files/fileGnssSignalBias.h"
#include "files/fileSatelliteModel.h"
#include "config/config.h"
#include "classes/eclipse/eclipse.h"
#include "classes/ephemerides/ephemerides.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "classes/parametrizationAcceleration/parametrizationAcceleration.h"
#include "misc/observation/variationalEquation.h"
#include "misc/observation/variationalEquationFromFile.h"
#include "gnss/gnss.h"
#include "gnss/gnssTransmitter.h"
#include "gnss/gnssParametrizationAntenna.h"

/**
* @defgroup gnssParametrizationTransmitterGroup GnssParametrizationTransmitter
* @brief Transmitter satellite constellations.
* @ingroup gnssGroup
* Parametrization of GNSS transmitter satellites.
* The interface is given by @ref GnssParametrizationTransmitter.
* An Instance can be created by @ref readConfig. */
/// @{

/***** TYPES ***********************************/

class GnssParametrizationTransmitter;
typedef std::shared_ptr<GnssParametrizationTransmitter> GnssParametrizationTransmitterPtr;

/***** CLASS ***********************************/

/** @brief Transmitter satellite constellations.
* Parametrization of GNSS transmitter satellites.
* An Instance of this class can be created by @ref readConfig. */
class GnssParametrizationTransmitter : public std::enable_shared_from_this<GnssParametrizationTransmitter>, public Gnss::Parametrization
{
public:

  // ======================

  class BiasModel
  {
  public:
    FileName                   fileNameOut;
    GnssType                   type;
    ParametrizationTemporalPtr temporal;
    Vector                     x;
    Vector                     bias; // for each epoch
    Gnss::ParameterIndex       indexParameter;
  };

  class TimeVariableBias
  {
  public:
    FileName inNameBias;
    GnssType type;
    Vector biases;
  };

  // ======================

  class Transmitter : public Gnss::Transmitter
  {
  protected:
    GnssParametrizationTransmitter *base;
    GnssStationInfo               transmitterInfo;
    GnssType                      type; // system + PRN
    Bool                          useAtAll;
    std::vector<Bool>             use;
    std::vector<Transform3d>      crf2arf;
    std::vector<TimeVariableBias> timeVariableBiases;

    // Parametrization clock error
    // ---------------------------
    Vector                            clock0, dclock; // [m]
    Double                            clockSigma;
    Matrix                            clockNoise;
    std::vector<Gnss::ParameterIndex> indexParameterClock;
    Gnss::ParameterIndex              indexParameterClockModel;

    // Parametrization of variational equations
    // ----------------------------------------
    VariationalEquationFromFile::ObservationEquation variationalEquation;   //!< Contains variational data: times, pos0, vel0, PosDesign, VelDesign, rotSat, rotEarth
    Gnss::ParameterIndex       indexParameterSatellite;                     //!< Index of unknown satellite related parameters in normal equation
    UInt                       countParameterSatellite;                     //!< Number of unknown satellite related parameters, includes stochastic pulses (x,y,z)
    Gnss::ParameterIndex       indexParameterSatelliteArc;                  //!< Index of unknown arc related parameters in normal equation
    UInt                       countParameterSatelliteArc;                  //!< Number of unknown arc related parameters, includes satellite state (6 parameter)
    std::vector<ParameterName> parameterNameSatellite;                      //!< Names of satellite related parameters
    std::vector<ParameterName> parameterNameSatelliteArc;                   //!< Names of arc related parameters
    Vector                     xOrbit;                                      //!< Estimated parameters (solar radiation pressure, stochastic pulses, initial state)
    SatelliteModelPtr          satelliteModel;                              //!< Satellite macro model from variational file

    // Parametrization signal bias
    // ---------------------------
    GnssSignalBias             signalBias;
    Gnss::ParameterIndex       indexParameterBias;
    Matrix                     Bias;  // Transformation parameterBiasType -> SignalBiasType
    std::vector<BiasModel>     biasModel;

  public:
    Transmitter() {}
   ~Transmitter() {}

    GnssType    PRN() const override {return type;} // prn + GnssType::SYSTEM
    std::string name() const override;
    Bool        useable    (UInt idEpoch=NULLINDEX) const override;
    void        disable    (UInt idEpoch=NULLINDEX);
    Vector3d    positionCoM(UInt idEpoch, const Time &time) const;           // Center of Mass in CRF
    Vector3d    position   (UInt idEpoch, const Time &time) const override;  // antenna reference point in CRF
    Vector3d    velocity   (UInt idEpoch, const Time &time) const override;  // in CRF
    Transform3d celestial2antennaFrame(UInt idEpoch, const Time &time) const override;  // CRF -> satellite left-handed antenna frame
    Double      clockError (UInt idEpoch, const Time &time) const override;  // error = transmitter clock time - system time
    std::vector<GnssType> definedTypes(UInt idEpoch) const override;         // transmitted signal types
    Vector      antennaVariations(UInt idEpoch, Angle azimut, Angle elevation, const std::vector<GnssType> &type) const override;
    Bool        isDesignMatrixTransmitter(const Gnss::NormalEquationInfo &normalEquationInfo, UInt idRecv, UInt idEpoch) const override;
    void        designMatrixTransmitter(const Gnss::NormalEquationInfo &normalEquationInfo, const Gnss::ObservationEquation &eqn, Gnss::DesignMatrix &A) const override;
    Bool        supportsIntegerAmbiguities(const Gnss::NormalEquationInfo &normalEquationInfo) const override;
    Bool        isCodeBiasEstimated       (const Gnss::NormalEquationInfo &normalEquationInfo) const override;
    Bool        isPhaseBiasEstimated      (const Gnss::NormalEquationInfo &normalEquationInfo) const override;

    friend class GnssParametrizationTransmitter;
  };

  // ======================

  // constructor
  GnssParametrizationTransmitter();
  virtual ~GnssParametrizationTransmitter() {}

  void init(const FileName &fileNameTransmitterList, const FileName &fileNameTransmitterInfo,
            const FileName &fileNameAntennaDef, const FileName &fileNameReceiverDef, UInt interpolationDegree);

  std::vector<Gnss::TransmitterPtr>     transmitters();
  std::vector<Gnss::ParametrizationPtr> parametrizations();

  virtual std::string name()   const = 0;
  virtual GnssType    system() const = 0;

  // parametrization
  // ---------------
  void   initIntervalTransmitter (Gnss::AnalysisType analysisType, const std::vector<Time> &times, const Time &timeMargin, Parallel::CommunicatorPtr comm) override;
  Bool   initIntervalDisabling   (Gnss::AnalysisType analysisType, const std::vector<Time> &times, const Time &timeMargin, Parallel::CommunicatorPtr comm) override;
  void   initIntervalLate        (Gnss::AnalysisType analysisType, const std::vector<Time> &times, const Time &timeMargin, Parallel::CommunicatorPtr comm) override;
  void   initParameter           (Gnss::NormalEquationInfo &normalEquationInfo) override;
  void   aprioriParameter        (const Gnss::NormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const override;
  void   observationEquationEpoch(const Gnss::NormalEquationInfo &normalEquationInfo, UInt idEpoch, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const override;
  void   observationEquation(const Gnss::NormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const override;
  Double updateParameter(const Gnss::NormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz, Bool printStatistics) override;
  void   writeResults(const Gnss::NormalEquationInfo &normalEquationInfo, const std::string &suffix) override;

  /** @brief creates an derived instance of this class. */
  static GnssParametrizationTransmitterPtr create(Config &config, const std::string &name);

protected:
  std::vector<std::shared_ptr<Transmitter>> transmitter;

  FileName fileNameTemplateInVariational;
  FileName fileNameTemplateInClock;
  FileName fileNameTemplateInSignal;

  FileName fileNameTemplateOutVariational;
  FileName fileNameTemplateOutOrbit;
  FileName fileNameTemplateOutParameter;
  FileName fileNameTemplateOutSignal;
  FileName fileNameTemplateOutClock;
  FileName fileNameOutTransmitterList;

  std::vector<Time> times;
  Bool              disableShadowEpochs;
  Bool              disablePostShadowRecoveryEpochs;
  std::vector<TimeVariableBias> timeVariableBiases;
  GnssAntennaDefinition::NoPatternFoundAction noPatternFoundAction;

  // position interpolation
  // ----------------------
  Polynomial polynomial;

  // Parametrization of clock errors
  // -------------------------------
  FileName                   fileNameOutClock;
  FileName                   fileNameClockNoise;
  enum class EstimateClockError {NONE, EPOCH, NOISEMODEL};
  Double                     sigmaClockZeroMean;
  EstimateClockError         estimateClockError;
  Bool                       estimateClockSigma;
  ParametrizationTemporalPtr clockModel;

  // Parametrization of variational equations
  // ----------------------------------------
  Bool                           estimateDynamicOrbit;
  Double                         minEstimableEpochsRatio;
  UInt                           integrationDegree;
  EphemeridesPtr                 ephemerides;
  ParametrizationAccelerationPtr parametrizationAcceleration;
  std::vector<Time>              stochasticPulse;
  Bool                           onlyEclipsingStochasticPulse;
  EclipsePtr                     eclipse;

  // Parametrization signal bias
  // ---------------------------
  Bool                   estimatePhaseBias;
  Bool                   estimateCodeBias;
  Bool                   integerAmbiguities;
  std::vector<BiasModel> biasModel;

  // Parametrization of antenna center variations
  // --------------------------------------------
  std::vector<GnssParametrizationAntennaPtr> antennaCenterVariations;

  friend class Transmitter;
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class GnssParametrizationTransmitter.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a var is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] var Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates Troposphere */
template<> Bool readConfig(Config &config, const std::string &name, GnssParametrizationTransmitterPtr &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

template<> Bool readConfig(Config &config, const std::string &name, GnssParametrizationTransmitter::BiasModel &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***********************************************/

#endif /* __GROOPS___ */
