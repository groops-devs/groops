/***********************************************/
/**
* @file gnssParametrizationReceiverLeo.h
*
* @brief GNSS for Low Earth Orbiter (LEO).
*
* @author Torsten Mayer-Guerr
* @author Norbert Zehentner
* @author Sebastian Strasser
* @date 2010-08-03
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPARAMETRIZATIONRECEIVERLEO__
#define __GROOPS_GNSSPARAMETRIZATIONRECEIVERLEO__

// Latex documentation
#ifdef DOCSTRING_GnssParametrizationReceiver
static const char *docstringGnssParametrizationReceiverLeo = R"(
\section{LowEarthOrbiter}\label{gnssParametrizationReceiverType:lowEarthOrbiter}

A single low-Earth orbiting (LEO) satellite with an onboard GNSS receiver. Orbit positions and receiver clock errors are estimated per epoch.

Attitude data must be provided via \configFile{inputfileStarCamera}{instrument}. If no attitude data is available from the satellite operator,
the star camera data can be simulated by using \program{SimulateStarCamera} or \program{SimulateStarCameraSentinel1}.

Antenna center variations can be estimated by defining \configClass{antennaCenter}{gnssParametrizationAntennaType}.

In a \config{preprocessing} step observation outliers are removed or downweighted, continuous tracks of phase observations are
defined for ambiguity parametrization, and cycle slips are detected.
)";
#endif

/***********************************************/

#include "files/fileGnssStationInfo.h"
#include "files/fileGnssSignalBias.h"
#include "gnss/gnss.h"
#include "gnss/gnssReceiver.h"
#include "gnss/gnssParametrizationAntenna.h"
#include "gnss/gnssParametrizationReceiver.h"

/***** CLASS ***********************************/

/** @brief GNSS for Low Earth Orbiter (LEO).
* @ingroup gnssParametrizationReceiverGroup
* @see GnssParametrizationReceiver */
class GnssParametrizationReceiverLeo : public std::enable_shared_from_this<GnssParametrizationReceiverLeo>,
                                       public GnssParametrizationReceiver, public Gnss::Receiver,  public Gnss::Parametrization
{
private:
  FileName fileNameOutOrbit, fileNameOutClock, fileNameOutCovariance, fileNameResiduals;
  FileName fileNameOutSignal;
  FileName fileNameObs, fileNameOrbit, fileNameStarCamera;
  FileName fileNameSignalBias;

  Time timeIntervalStart, timeIntervalEnd;

  GnssStationInfo stationInfo;  // antenna reference point, ...
  GnssSignalBias  signalBias;   // signal Bias (DCBs, ...)

  // values at epochs
  // ----------------
  std::vector<Vector3d>    pos0;         // aproximate value (reference point in CRF, e.g. CoM)
  std::vector<Vector3d>    displacement; // antenna center in CRF relative to reference point
  std::vector<Vector3d>    dpos;         // estimated correction relative to pos0
  std::vector<Tensor3d>    cov;          // epoch covariances
  std::vector<Vector3d>    vel;          // velocity
  std::vector<Transform3d> crf2sat;      // transformation CRF -> satellite frame

  // Preprocessing
  // -------------
  std::vector<GnssType> useType, ignoreType;
  Angle                 cutOff, elevationTrackMinimum;
  Double                codeMaxPosDiff;
  UInt                  minObsCountPerTrack;
  Double                windowROTI;
  mutable VariableList  varList;
  ExpressionVariablePtr exprSigmaPhase;
  ExpressionVariablePtr exprSigmaCode;
  GnssAntennaDefinition::NoPatternFoundAction noPatternFoundAction;

  // Parametrization
  // ---------------
  std::vector<Gnss::ParameterIndex>     indexParameterEpoch;
  Gnss::ParameterIndex                  indexParameterBias;
  Bool                                  estimatePhaseBias;
  Bool                                  estimateCodeBias;
  Bool                                  integerAmbiguities;
  Double                                wavelengthFactor_;
  Matrix                                Bias; // Transformation parameterBiasType -> SignalBiasType
  std::vector<GnssParametrizationAntennaPtr> antennaCenterVariations;

public:
  GnssParametrizationReceiverLeo(Config &config);
  virtual ~GnssParametrizationReceiverLeo() {}

  std::vector<Gnss::ParametrizationPtr> additionalParametrizations() const {return std::vector<Gnss::ParametrizationPtr>(antennaCenterVariations.begin(), antennaCenterVariations.end());}

  // GnssReceivers
  // -------------
  std::vector<Gnss::ReceiverPtr>        receivers()        override {return std::vector<Gnss::ReceiverPtr>({shared_from_this()});}
  std::vector<Gnss::ParametrizationPtr> parametrizations() override;

  // Parametrization
  // ---------------
  void   initIntervalReceiver (Gnss::AnalysisType analysisType, const std::vector<Time> &times, const Time &timeMargin, Parallel::CommunicatorPtr comm) override;
  Bool   initIntervalDisabling(Gnss::AnalysisType analysisType, const std::vector<Time> &times, const Time &timeMargin, Parallel::CommunicatorPtr comm) override;
  void   initIntervalLate     (Gnss::AnalysisType analysisType, const std::vector<Time> &times, const Time &timeMargin, Parallel::CommunicatorPtr comm) override;
  void   initParameter(Gnss::NormalEquationInfo &normalEquationInfo) override;
  void   aprioriParameter(const Gnss::NormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const override;
  Double updateParameter(const Gnss::NormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz, Bool printStatistics) override;
  void   updateCovariance(const Gnss::NormalEquationInfo &normalEquationInfo, const MatrixDistributed &covariance) override;
  void   writeResults(const Gnss::NormalEquationInfo &normalEquationInfo, const std::string &suffix) override;

  // receiver
  // --------
  std::string name() const  override {return stationInfo.markerName;}
  Vector3d    position         (UInt idEpoch) const override {return pos0.at(idEpoch) + displacement.at(idEpoch) + dpos.at(idEpoch);}
  Vector3d    velocity         (UInt idEpoch) const override {return vel.at(idEpoch);}
  Transform3d celestial2localFrame(UInt idEpoch) const override {return crf2sat.at(idEpoch);}
  Transform3d local2antennaFrame  (UInt idEpoch) const override;
  Double      troposphere      (UInt /*idEpoch*/, Angle /*azimut*/, Angle /*elevation*/) const override {return 0.0;}
  std::vector<GnssType> definedTypes(UInt idEpoch) const override;
  Vector      antennaVariations(UInt idEpoch, Angle azimut, Angle elevation, const std::vector<GnssType> &type) const override;
  Vector      accuracy(UInt idEpoch, Angle azimut, Angle elevation, const std::vector<GnssType> &type) const override;
  Angle       elevationCutOff() const override {return cutOff;}
  Bool        isEpochEstimable(Gnss::AnalysisType analysisType, UInt idEpoch) const override;
  Bool        isDesignMatrixReceiver(const Gnss::NormalEquationInfo &normalEquationInfo, UInt idTrans, UInt idEpoch) const override;
  void        designMatrixReceiver(const Gnss::NormalEquationInfo &normalEquationInfo, const Gnss::ObservationEquation &eqn, Gnss::DesignMatrix &A) const override;
  Bool        supportsIntegerAmbiguities(const Gnss::NormalEquationInfo &/*normalEquationInfo*/) const override {return integerAmbiguities;}
  Double      wavelengthFactor() const override {return wavelengthFactor_;}
  Bool        isCodeBiasEstimated (const Gnss::NormalEquationInfo &normalEquationInfo) const override;
  Bool        isPhaseBiasEstimated(const Gnss::NormalEquationInfo &normalEquationInfo) const override;
};

/***********************************************/

#endif /* __GROOPS___ */
