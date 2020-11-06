/***********************************************/
/**
* @file gnssParametrizationReceiverStationNetwork.h
*
* @brief GNSS ground station network.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2010-08-19
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPARAMETRIZATIONRECEIVERSTATIONNETWORK__
#define __GROOPS_GNSSPARAMETRIZATIONRECEIVERSTATIONNETWORK__

// Latex documentation
#ifdef DOCSTRING_GnssParametrizationReceiver
static const char *docstringGnssParametrizationReceiverStationNetwork = R"(
\section{StationNetwork}\label{gnssParametrizationReceiverType:stationNetwork}

A network of GNSS ground stations is defined via \configFile{inputfileStationList}{stringList}.
All input and output files except \configFile{outputfileGriddedPosition}{griddedData},
\configFile{outputfileUsedStationList}{stringList}, \configFile{inputfileAntennaDefinition}{gnssAntennaDefinition},
\configFile{inputfileReceiverDefinition}{gnssReceiverDefinition}, and
\configFile{inputfileAccuracyDefinition}{gnssAntennaDefinition} are read/written for each station.
The file name is interpreted as a template with the variable \verb|{station}| being replaced by the station name.

If \config{kinematicPositionEstimation}=\verb|yes| station positions are estimated every epoch,
otherwise they are parameterized according to \configClass{positionEstimation}{parametrizationTemporalType}.

In case you estimate the station network together with a GNSS constellations, \config{noNetTranslationSigma} and \config{noNetRotationSigma}
can be set to remove the rank deficiency. No-net constraints are only applied to stations defined in \configFile{noNetStationList}{stringList} that have
precise station coordinates provided via \configFile{inputfileStationPosition}{instrument}. Those will be used instead of the approximate
coordinates from \configFile{inputfileStationInfo}{gnssStationInfo}. In case you want to align to an ITRF/IGS reference frame, precise coordinates can be
generated via combining \program{Sinex2StationPosition}, \program{Sinex2StationDiscontinuities}, and \program{Sinex2StationPostSeismicDeformation}.

A priori tropospheric correction is handled by a \configClass{troposphere}{troposphereType} model (e.g. Vienna Mapping Functions 3).
Additional parameters for zenith wet delay and gradients can be set up via \configClass{troposphereWetEstimation}{parametrizationTemporalType}
and \configClass{troposphereGradientEstimation}{parametrizationTemporalType}. These parameters can be soft-constrained using
\configClass{constraints}{gnssParametrizationConstraintsType} in \program{GnssProcessing}.

Antenna center variations can be estimated by defining \configClass{antennaCenter}{gnssParametrizationAntennaType}.
It is possible to estimate patterns for each individual station or grouped by antenna model.
Antenna center parameters can be constrained by adding \configClass{constraints}{gnssParametrizationConstraintsType} in \program{GnssProcessing}.

The effects of loading and tidal deformation on station positions can be corrected for via \configClass{deformation}{gravityfieldType} and
\configClass{tides}{tidesType}, respectively. These typically include:
\begin{itemize}
  \item \configClass{tides:earthTide}{tidesType:earthTide}: Earth tidal deformations (IERS conventions)
  \item \configClass{tides:doodsonHarmonicTide}{tidesType:doodsonHarmonicTide}: ocean tidal deformations
        (e.g. FES 2014b, \config{minDegree}=\verb|1|)
  \item \configClass{tides:doodsonHarmonicTide}{tidesType:doodsonHarmonicTide}: atmospheric tidal deformation
        (e.g. AOD1B RL06, \config{minDegree}=\verb|1|)
  \item \configClass{tides:poleTide}{tidesType:poleTide}: pole tidal deformations (IERS conventions)
  \item \configClass{tides:poleOceanTide}{tidesType:oceanPoleTide}: ocean pole tidal deformations (IERS conventions)
\end{itemize}

Each station goes through a \config{preprocessing} step individually, where observation outliers are removed or downweighted,
continuous tracks of phase observations are defined for ambiguity parametrization, cycle slips are detected, and stations are
disabled if they do not fulfill certain requirements.
)";
#endif

/***********************************************/

#include "config/config.h"
#include "files/fileGnssSignalBias.h"
#include "files/fileGnssStationInfo.h"
#include "classes/ephemerides/ephemerides.h"
#include "classes/troposphere/troposphere.h"
#include "classes/gravityfield/gravityfield.h"
#include "classes/tides/tides.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "gnss/gnss.h"
#include "gnss/gnssReceiver.h"
#include "gnss/gnssParametrizationAntenna.h"
#include "gnss/gnssParametrizationEarthRotation.h"
#include "gnss/gnssParametrizationReceiver.h"

/***** CLASS ***********************************/

/** @brief GNSS for permanent stations.
* @ingroup gnssParametrizationReceiverGroup
* @see GnssParametrizationReceiver */
class GnssParametrizationReceiverStationNetwork : public std::enable_shared_from_this<GnssParametrizationReceiverStationNetwork>,
                                           public GnssParametrizationReceiver, public Gnss::Parametrization
{
public:

  class Receiver : public Gnss::Receiver
  {
  public:
    GnssParametrizationReceiverStationNetwork  *base;
    std::string           stationName;
    GnssStationInfo       stationInfo;
    GnssSignalBias        signalBias;
    UInt                  countAlternatives; // number of following alternative stations
    UInt                  idTropo;

    std::vector<Double>   zenitDelayWet;
    std::vector<Double>   gradientX, gradientY;
    Transform3d           lnof2trf; // north, east, up -> TRF

    // positions in TRF
    Vector3d              approxPosition;  // approximate value
    std::vector<Vector3d> posDisplacement; // tides and loading variations + antenna center (ARP)
    std::vector<Vector3d> dpos;            // estimated refinements

    // parametrization
    // ---------------
    Vector               xPos;
    Vector               xTropoWet;
    Vector               xTropoGradient;
    std::vector<Gnss::ParameterIndex> indexParameterEpoch;
    Gnss::ParameterIndex indexParameterPosition;
    Gnss::ParameterIndex indexParameterTropoWet;
    Gnss::ParameterIndex indexParameterTropoGradient;
    Gnss::ParameterIndex indexParameterBias;
    Matrix               Bias; // Transformation parameterBiasType -> SignalBiasType

    static void designMatrixTemporal(ParametrizationTemporalPtr parametrization, const Time &time, const_MatrixSliceRef B,
                                     const Gnss::ParameterIndex &index, Gnss::DesignMatrix &A);
  public:
    Receiver() {}
   ~Receiver() {}

    std::string name() const override;
    Vector3d    position                  (UInt idEpoch) const override;
    Vector3d    velocity                  (UInt idEpoch) const override;
    Transform3d celestial2localFrame      (UInt idEpoch) const override;
    Transform3d local2antennaFrame        (UInt idEpoch) const override;
    Double      troposphere               (UInt idEpoch, Angle azimut, Angle elevation) const override;
    std::vector<GnssType> definedTypes    (UInt idEpoch) const override;
    Vector      antennaVariations         (UInt idEpoch, Angle azimut, Angle elevation, const std::vector<GnssType> &type) const override;
    Vector      accuracy                  (UInt idEpoch, Angle azimut, Angle elevation, const std::vector<GnssType> &type) const override;
    Angle       elevationCutOff() const override;
    Bool        isEpochEstimable          (Gnss::AnalysisType analysisType, UInt idEpoch) const override;
    Bool        isDesignMatrixReceiver    (const Gnss::NormalEquationInfo &normalEquationInfo, UInt idTrans, UInt idEpoch) const override;
    void        designMatrixReceiver      (const Gnss::NormalEquationInfo &normalEquationInfo, const Gnss::ObservationEquation &eqn, Gnss::DesignMatrix &A) const override;
    Bool        supportsIntegerAmbiguities(const Gnss::NormalEquationInfo &normalEquationInfo) const override;
    Bool        isCodeBiasEstimated       (const Gnss::NormalEquationInfo &normalEquationInfo) const override;
    Bool        isPhaseBiasEstimated      (const Gnss::NormalEquationInfo &normalEquationInfo) const override;

    friend class GnssParametrizationReceiverStationNetwork;
  };

  // constructor
  GnssParametrizationReceiverStationNetwork(Config &config);
  virtual ~GnssParametrizationReceiverStationNetwork() {}

  // GnssReceivers
  // -------------
  std::vector<Gnss::ReceiverPtr>        receivers()        override;
  std::vector<Gnss::ParametrizationPtr> parametrizations() override;

  void   initIntervalReceiver    (Gnss::AnalysisType analysisType, const std::vector<Time> &times, const Time &timeMargin, Parallel::CommunicatorPtr comm) override;
  Bool   initIntervalDisabling   (Gnss::AnalysisType analysisType, const std::vector<Time> &times, const Time &timeMargin, Parallel::CommunicatorPtr comm) override;
  void   initIntervalLate        (Gnss::AnalysisType analysisType, const std::vector<Time> &times, const Time &timeMargin, Parallel::CommunicatorPtr comm) override;
  void   initParameter           (Gnss::NormalEquationInfo &normalEquationInfo) override;
  void   aprioriParameter        (const Gnss::NormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const override;
  void   observationEquation     (const Gnss::NormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const override;
  Double updateParameter         (const Gnss::NormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz, Bool printStatistics) override;
  void   writeResults            (const Gnss::NormalEquationInfo &normalEquationInfo, const std::string &suffix) override;

private:
  std::vector<std::shared_ptr<Receiver>> receiver;

  FileName              fileNameGrid;
  FileName              fileNamePosition;
  FileName              fileNamePositionSeries;
  FileName              fileNameResiduals;
  FileName              fileNameOutSignal;
  FileName              fileNameUsedStationList;
  FileName              fileNameOutClock;
  FileName              fileNameOutTropo;

  FileName              fileNameStationPosition;
  FileName              fileNameObs;
  FileName              fileNameSignalBias;
  std::vector<Time>     times;

  FileName              fileNameNoNetStationList;
  Double                sigmaNoNetTranslation, sigmaNoNetRotation;
  Vector                netTranslation, netRotation;
  Vector                noNetStations;

  // Preprocessing
  // -------------
  FileName              fileNameTrackBefore, fileNameTrackAfter;
  UInt                  maxStationCount;
  std::vector<GnssType> useType, ignoreType;
  Angle                 elevationCutOff, elevationTrackMinimum;
  Double                codeMaxPosDiff;
  UInt                  minObsCountPerTrack;
  Double                minEstimableEpochsRatio;
  Double                denoisingLambda;
  UInt                  tecWindowSize;
  Double                tecSigmaFactor;
  GnssAntennaDefinition::NoPatternFoundAction noPatternFoundAction;

  // Parametrization of positions
  // ----------------------------
  Bool                       isKinematicPosition;
  ParametrizationTemporalPtr parametrizationPosition;

  // Parametrization signal bias
  // ---------------------------
  Bool estimatePhaseBias;
  Bool estimateCodeBias;
  Bool integerAmbiguities;

  // Tropospheric delay
  // ------------------
  TropospherePtr             troposphere;
  ParametrizationTemporalPtr parametrizationTropoWet;
  ParametrizationTemporalPtr parametrizationTropoGradient;

  // Parametrization of antenna center variations
  // --------------------------------------------
  std::vector<GnssParametrizationAntennaPtr> antennaCenterVariations;

  // tides & loading deformation
  // ---------------------------
  Vector          hn,ln;
  EphemeridesPtr  ephemerides;
  GravityfieldPtr gravityfield;
  TidesPtr        tides;
};

/***********************************************/

#endif /* __GROOPS___ */
