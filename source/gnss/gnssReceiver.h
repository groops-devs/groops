/***********************************************/
/**
* @file gnssReceiver.h
*
* @brief GNSS receiver.
*
* @author Torsten Mayer-Guerr
* @author Norbert Zehentner
* @author Sebastian Strasser
* @date 2013-06-28
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSRECEIVER__
#define __GROOPS_GNSSRECEIVER__

#include "gnss/gnss.h"
#include "files/fileInstrument.h"

/** @addtogroup gnssGroup */
/// @{

/***** CLASS ***********************************/

/** @brief Abstract class for GNSS receiver.
* eg. permanent stations or LEOs. */
class Gnss::Receiver
{
  Gnss   *_gnss;    // is set by registerGnss
  UInt    _idRecv;  // is set by registerGnss

  Bool                  useAtAll;
  std::vector<Bool>     use;
  std::vector<Double>   clk;    // clock error
  std::vector<Time>     times;  // observed time

  std::vector<GnssType> sigma0Type;
  std::vector<Double>   sigma0Factor;

  Double                observationSampling;
  std::vector<Gnss::Observation> obsMem;
  std::vector<std::vector<Observation*>> observations_; //!< observations at receiver (for each epoch, for each transmitter)

protected:
  Double huber;
  Double huberPower;

public:
  // Types
  // -----
  enum WeightingType {NONE, INDIVIDUAL, GROUPED, GROUPEDPHASES};

  // public variables
  // ----------------
  std::vector<Track*> track; //!< tracking phase observations

// ========================================================

/// Constructor.
Receiver();

/// Destructor.
virtual ~Receiver();

/** @brief delete observations and tracks. */
void free();

/** @brief is called by Gnss::init(). */
void registerGnss(Gnss *gnss, UInt idRecv) {_gnss = gnss; _idRecv = idRecv;}

/** @brief Reference to the complete GNSS system. */
Gnss &gnss()   const;

/** @brief Identify number in the GNSS system. */
UInt  idRecv() const {return _idRecv;}

/** @brief Name of the station or satellite */
virtual std::string name() const = 0;

// ========================================================

/** @brief Initialize useable() and clockError().
* Must be called at one node only. */
void initInterval(const std::vector<Time> &times);

/** @brief Is the receiver usable at given epoch (or all epochs).
* Can be useable at one node only. */
Bool useable(UInt idEpoch=NULLINDEX) const {return useAtAll && ((idEpoch == NULLINDEX) || use.at(idEpoch));}

/** @brief Set the receiver usable at given epoch (or all epochs). */
void disable(UInt idEpoch=NULLINDEX);

// ========================================================

/** @brief Clock error.
* error = receiver clock time - system time [s] */
Double clockError(UInt idEpoch) const {return clk.at(idEpoch);}

/** @brief set clock error.
* error = receiver clock time - system time [s] */
virtual void updateClockError(UInt idEpoch, Double deltaClock) {clk.at(idEpoch) += deltaClock;}

/** @brief Receiver time of received signals.
* (receiver clock time disturbed by clock error) */
Time timeReceiver(UInt idEpoch) const {return times.at(idEpoch);}

/** @brief Estimated system time of received signals.
* (receiver clock time - clock error) */
Time timeCorrected(UInt idEpoch) const {return times.at(idEpoch) - seconds2time(clockError(idEpoch));}

// ========================================================

/** @brief antenna reference point in CRF. */
virtual Vector3d position(UInt idEpoch) const = 0;

/** @brief velocity in CRF [m/s]. */
virtual Vector3d velocity(UInt idEpoch) const = 0;

/** @brief Rotation from celestial reference frame (CRF) to local/body frame (north, east, up or vehicle system). */
virtual Transform3d celestial2localFrame(UInt idEpoch) const = 0;

/** @brief Rotation from local/body frame to left-handed antenna system (usually north, east, up). */
virtual Transform3d local2antennaFrame(UInt idEpoch) const = 0;

/** @brief elevation cut off angle. */
virtual Angle elevationCutOff() const = 0;

/** @brief Tropspheric delay. */
virtual Double troposphere(UInt idEpoch, Angle azimut, Angle elevation) const = 0;

/** @brief Observed signal types. Empty if no GNSS receiver definition was provided. */
virtual std::vector<GnssType> definedTypes(UInt idEpoch) const = 0;

/** @brief Direction dependent corrections.
* Antenna center variations + signal biases.
* observed range = range (ARPs of transmitter and receiver) + antennaVariations.
* @a azmiut and @a elevation must be given in the antenna frame (left-handed). */
virtual Vector antennaVariations(UInt idEpoch, Angle azimut, Angle elevation, const std::vector<GnssType> &type) const = 0;

/** @brief Direction (and other parameters) dependent standard deviation.
* @a azmiut and @a elevation must be given in the antenna frame (left-handed). */
virtual Vector accuracy(UInt idEpoch, Angle azimut, Angle elevation, const std::vector<GnssType> &type) const = 0;

/** @brief Transformation matrix for observed (composed) types from orignal transmitted types.
* E.g. C2DG = C1CG - C1WG + C2WG.
* Returns the @a typesTrans and the transformation matrix @a A (dimension: types.size() times typesTrans.size()). */
virtual void signalComposition(UInt /*idEpoch*/, const std::vector<GnssType> &types, std::vector<GnssType> &typesTrans, Matrix &A) const {Gnss::defaultSignalComposition(types, typesTrans, A);}

// ========================================================

/** @brief Enough observations? */
virtual Bool isEpochEstimable(AnalysisType analysisType, UInt idEpoch) const = 0;

/** @brief Enough estimable epochs? @p minEstimableEpochsRatio must be between 0 and 1. Factors in receiver observation sampling. */
Bool isReceiverEstimable(Double minEstimableEpochsRatio) const;

/** @brief Parameters of receiver of this observation. */
virtual Bool isDesignMatrixReceiver(const NormalEquationInfo &normalEquationInfo, UInt idTrans, UInt idEpoch) const = 0;

/** @brief Fill in the design matrix with parameters of receiver. */
virtual void designMatrixReceiver(const NormalEquationInfo &normalEquationInfo, const ObservationEquation &eqn, DesignMatrix &A) const = 0;

/** @brief Factor to account for half-wavelength observations (collected by codeless squaring techniques).  */
virtual Double wavelengthFactor() const {return 1.;}

/** @brief Receiver is able to track full cycle integer ambiguities. */
virtual Bool supportsIntegerAmbiguities(const NormalEquationInfo &normalEquationInfo) const = 0;

/** @brief Are code biases estimated?. */
virtual Bool isCodeBiasEstimated(const NormalEquationInfo &normalEquationInfo) const = 0;

/** @brief Are (float) phase biases estimated?. */
virtual Bool isPhaseBiasEstimated(const NormalEquationInfo &normalEquationInfo) const = 0;

// ========================================================

/** @brief All observations between receiver and transmitter at one epoch. */
Observation *observation(UInt idTrans, UInt idEpoch) const;

/** @brief Max. observed epoch id +1. */
UInt idEpochSize() const {return observations_.size();}

/** @brief Max. observed transmitter id+1 at @a idEpoch. */
UInt idTransmitterSize(UInt idEpoch) const {return observations_.at(idEpoch).size();}

/** @brief Counts observations between receiver and transmitter(s). */
UInt countObservations(UInt idTrans, UInt idEpochStart, UInt idEpochEnd) const;

/** @brief Number and types of observations. */
void countObservationsPerType(UInt idTrans, UInt idEpochStart, UInt idEpochEnd, std::vector<GnssType> &type, std::vector<UInt> &count) const;

/** @brief Types of observations. */
std::vector<GnssType> observationsTypes(UInt idTrans) const;

/** @brief Computes residual statistics. */
void residualsStatistics(UInt idTrans, const std::vector<UInt> &idEpochs,
                         std::vector<GnssType> &type, std::vector<Double> &ePe, std::vector<Double> &redundancy,
                         std::vector<UInt> &obsCount, std::vector<UInt> &outlierCount) const;

/** @brief Computes weights from residuals. */
void computeWeightsFromResiduals(const std::vector<UInt> &idEpochs, WeightingType computeWeights, Gnss::Receiver::WeightingType adjustSigma0);

/** @brief Returns sigma factor for observation type. */
Double sigmaFactor(GnssType type) const;

/** @brief Returns residuals as satellite arc. */
GnssReceiverArc residuals(const std::vector<UInt> &idEpochs) const;

// Preprocessing
// -------------
  /** @brief Compute observation equations (reduced observations). */
  class ObservationEquationList
  {
    std::vector<std::vector<std::unique_ptr<ObservationEquation>>> eqn;

  public:
    ObservationEquationList();
    ObservationEquationList(const Gnss &gnss, const Receiver &receiver, AnalysisType analysisType);

    ObservationEquation *operator()(UInt idTrans, UInt idEpoch) const;
    void deleteObservationEquation(UInt idTrans, UInt idEpoch);
  };

  friend class ObservationEquationList;

protected:
/** @brief Reads observations from a file. Member variable @a times must be set. */
void readObservations(InstrumentFile &fileReceiver, const std::vector<GnssType> &useType, const std::vector<GnssType> &ignoreType, const Time &timeMargin);

/** @brief Initializes observations. Receiver and Transmitter positions, orientations, ... must be initialized beforehand. */
void initObservation(Gnss::AnalysisType analysisType);

/** @brief Delete observation. */
void deleteObservation(UInt idTrans, UInt idEpoch);

/** @brief Delete observations that don't match the types from receiver definition. Does nothing if no receiver definition is provided. */
Bool deleteUndefinedObservations();

/** @brief Set observations of an epoch to unusable if there are not enough observations to estimate a solution in the epoch. */
Bool deleteObservationsOfInestimableEpochs(Gnss::AnalysisType analysisType);

/** @brief Estimate coarse receiver clock errors from a Precise Point Positioning (PPP) code solution.
* If @p estimateKinematicPosition is TRUE, the receiver position is estimated at each epoch, otherwise it is estimated once for all epochs.*/
void estimateInitialClockErrorFromCodeObservations(Double maxPosDiff, Bool estimateKinematicPosition=TRUE);

/** @brief Disable epochs if reduced code observations @p eqn exceed @p threshold (e.g. 100+ km for code cycle slips). */
void disableEpochsWithGrossCodeObservationOutliers(ObservationEquationList &eqn, Double threshold);

/** @brief Create tracks with continously identical phase observations.
* Tracks may contain gaps but must contain observations in at least @p minObsCountPerTrack epochs.
* Extra types are included (e.g. L5*G), but tracks must have at least two others phases at different frequencies.*/
void createTracks(UInt minObsCountPerTrack, const std::vector<GnssType> &extraTypes={});

/** @brief Write TEC and Melbourne-Wuebbena like combinations for tracks to file @p fileName. One file per track (file name variables: prn, trackTimeStart, trackTimeEnd). */
void writeTracks(const FileName &fileName, VariableList varList, const ObservationEquationList &eqnList);

/** @brief delete track and all related obersvations. */
void deleteTrack(ObservationEquationList &eqn, UInt idTrack);

/** @brief Split a @p track at @p idEpochSplit into two new tracks. Shortens the old track and returns the new track. Updates observation and observation equation @p eqn track assignments. */
Track *splitTrack(ObservationEquationList &eqn, Gnss::Track *track, UInt idEpochSplit);

/** @brief Removes tracks that never exceed @p minElevation (in radian). */
void removeLowElevationTracks(ObservationEquationList &eqn, Angle minElevation);

/** @brief Iteratively removes epochs with too few observed transmitters and subsequently removes tracks that don't have @p minObsCountPerTrack observations anymore. */
Bool removeTracksWithInsufficientObservationCount(ObservationEquationList &eqn, UInt minObsCountPerTrack, Gnss::AnalysisType analysisType);

/** @brief Track outlier detection based on robust least squares estimation. Downweights outliers and prereduces observations by estimated integer ambiguities. */
void trackOutlierDetection(const ObservationEquationList &eqn, const std::vector<GnssType> &ignoreTypes={});

/** @brief Robust least squares estimation for a single epoch. */
static Matrix robustLeastSquares(Double huber, Double power, const_MatrixSliceRef A, const_MatrixSliceRef l,
                                 const Vector &sigma0, Vector &sigma);

/** @brief Robust least squares estimation for multiple epochs in a single system of equations. */
static Matrix robustLeastSquares(Double huber, Double power,
                                 const std::vector<Matrix> &A, const std::vector<Matrix> &l, const std::vector<UInt> &obsCount,
                                 const Vector &sigma0, Vector &sigma);

/** @brief Splits tracks at detected cycle slips.
* Based on all Melbourne-Wuebbena like combinations.
* Solves the total variation regularized least-squares problem (total variation denoising).
* In a second step smoothness of TEC is evaluated using a moving window peak/outlier detection based on AR model residuals.
* Gaps in a track causes no problems.
* ASSUMPTION: All epochs contain the same phase observations.
* ATTENTION: Different number of code observations can generate apparent cycle slips as code biases are ignored.
* ATTENTION: GPS L5 observations are handled separately due to temporal changing bias.
* @param eqnList Observation equations (reduced observations).
* @param lambda Regularization parameter (@see @a totalVariationDenoising) (e.g. @p lambda = 5 for GPS ground stations).
* @param windowSize Size of the moving window used for the TEC smoothness evaluation. If 0, TEC is not analyzed.
* @param tecSigmaFactor Factor applied to moving standard deviation of AR model residuals to determine threshold for peak/outlier detection. */
void cycleSlipsDetection(ObservationEquationList &eqnList, Double lambda, UInt windowSize, Double tecSigmaFactor);
void cycleSlipsDetection(ObservationEquationList &eqnList, Track *track, Double lambda, UInt windowSize, Double tecSigmaFactor);

/** @brief repair cycle slip differences at same frequencies (e.g. between C1CG and C1WG).
* Allows to reduce the number of integer ambiguities.
* ASSUMPTION: All epochs of each track contain the same phase observations.
* @param eqnList Observation equations (reduced observations). */
void cycleSlipsRepairAtSameFrequency(ObservationEquationList &eqnList);


/** @brief Compute ROTI (Rate of TEC index) each track and add as peuso observations.
* GnssType R**+PRN, e.g. R**G32. TEC is estimated from all phase observations, the rate dTEC/dt as single difference,
* ROTI is the standard deviation at each time step of of all dTEC/dt in @a windowSize [t-windowSize/2, t+windowSize/2).
* ASSUMPTION: All epochs contain the same phase observations.
* @param eqnList Observation equations (reduced observations).
* @param windowSize [seconds] Size of the moving window. */
void addRotiPseudoObservations(ObservationEquationList &eqnList, Double windowSize);

public:
/** @brief Total variation denoising.
* Solves the total variation regularized least-squares problem.
* Laurent Condat. A Direct Algorithm for 1D Total Variation Denoising. IEEE Signal Processing Letters,
* Institute of Electrical and Electronics Engineers, 2013, 20 (11), pp.1054-1057. DOI: 10.1109/LSP.2013.2278339.
* @param y noisy input time series.
* @param lambda Regularization parameter.
* @return denoised time series. */
static std::vector<Double> totalVariationDenoising(const std::vector<Double> &y, Double lambda);
}; //class Gnss::Receiver

/***** CLASS ***********************************/

/** @brief Abstract class for ambiguities. */
class Gnss::Ambiguity
{
public:
  virtual ~Ambiguity() {}
  virtual Vector ambiguities(const std::vector<GnssType> &type) const = 0;
};

/***** CLASS ***********************************/

/** @brief Class for tracking phases along epochs. */
class Gnss::Track
{
public:
  ReceiverPtr           receiver;
  TransmitterPtr        transmitter;
  UInt                  idEpochStart, idEpochEnd;
  std::vector<GnssType> types;
  std::shared_ptr<Ambiguity> ambiguity;

  Track(ReceiverPtr _receiver, TransmitterPtr _transmitter, UInt _idEpochStart, UInt _idEpochEnd, const std::vector<GnssType> &_types) :
    receiver(_receiver), transmitter(_transmitter), idEpochStart(_idEpochStart), idEpochEnd(_idEpochEnd), types(_types) {}

  UInt countObservations() const;
  void removeAmbiguitiesFromObservations(const std::vector<GnssType> &type, const std::vector<Double> &value);

  /** @brief Form Melbourne-Wuebbena like linear combinations of phase observations.
  * ASSUMPTION: All epochs contain the same phase observations.
  * ATTENTION: Different number of code observations can generate apparent cycle slips as code biases are ignored.
  * ATTENTION: GPS L5 observations are ignored due to temporally changing bias.
  * @param eqnList Observation equations (reduced observations).
  * @param[out] idEpochs Epoch IDs for linear combinations.
  * @param[out] estimates Linear combinations per epoch, ordered from most inaccurate to most accurate.
  * @param[out] typesPhase Phase observation types.
  * @param[out] Bias Transformation matrix from phase decorrelation based on @p typePhase. */
  void combinations(const Gnss::Receiver::ObservationEquationList &eqnList, std::vector<UInt> &idEpochs,
                    std::vector<Vector> &estimates, std::vector<GnssType> &typesPhase, Matrix &Bias) const;

  /** @brief Determine range and TEC based on phase observations only.
  * ASSUMPTION: All epochs contain the same phase observations.
  * ATTENTION: GPS L5 observations are ignored due to temporally changing bias.
  * @param eqnList Observation equations (reduced observations).
  * @param typesPhase Phase observation types.
  * @param idEpochs Epoch IDs for linear combinations.
  * @param[out] range Range estimates at @p idEpochs in meters.
  * @param[out] tec TEC estimates at @p idEpochs in TEC units. */
  void rangeAndTecFromPhase(const Gnss::Receiver::ObservationEquationList &eqnList, const std::vector<GnssType> &typesPhase,
                            const std::vector<UInt> &idEpochs, std::vector<Double> &range, std::vector<Double> &tec) const;
};

/***********************************************/

/// @}

#endif /* __GROOPS___ */
