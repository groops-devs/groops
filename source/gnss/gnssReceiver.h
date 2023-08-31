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

#include "files/fileInstrument.h"
#include "classes/noiseGenerator/noiseGenerator.h"
#include "gnss/gnssObservation.h"
#include "gnss/gnssTransceiver.h"
#include "gnss/gnssTransmitter.h"

/** @addtogroup gnssGroup */
/// @{

/***** TYPES ***********************************/

class GnssReceiver;
class GnssTrack;
typedef std::shared_ptr<GnssReceiver> GnssReceiverPtr;
typedef std::shared_ptr<GnssTrack>    GnssTrackPtr;

/***** CLASS ***********************************/

/** @brief Abstract class for GNSS receiver.
* eg. permanent stations or LEOs. */
class GnssReceiver : public GnssTransceiver
{
  std::vector<GnssObservation> obsMem;
  std::vector<std::vector<GnssObservation*>> observations_; // observations at receiver (for each epoch, for each transmitter)

  void copyObservations2ContinuousMemoryBlock();

public:
  // public variables
  // ----------------
  Bool                      isMyRank_, isEarthFixed_;
  std::vector<Time>         times;
  std::vector<Double>       clk;
  std::vector<Vector3d>     pos, vel;  // regularized marker pos in global system
  std::vector<Vector3d>     offset;    // pos to ARF in local system
  std::vector<Transform3d>  global2local, global2antenna;
  std::vector<GnssTrackPtr> tracks;   // tracking phase observations
  Double                    observationSampling;
  Bool                      integerAmbiguities;
  Double                    wavelengthFactor; // Factor to account for half-wavelength observations (collected by codeless squaring techniques).
  std::vector<std::string>  preprocessingInfos;
  std::string               disableReason;

  GnssReceiver(Bool isMyRank, Bool isEarthFixed, const Platform &platform,
               GnssAntennaDefinition::NoPatternFoundAction noPatternFoundAction, const Vector &useableEpochs,
               Bool integerAmbiguities, Double wavelengthFactor);

  /// Destructor.
 ~GnssReceiver() {}

  /** @brief Identify number in the GNSS system. */
  UInt  idRecv() const {return id_;}

  /** @brief Disable given epoch (or all epochs). */
  void disable(UInt idEpoch, const std::string &reason) override;

  /** @brief Disable receiver completely. */
  void disable(const std::string &reason) override;

  Bool isMyRank() const {return isMyRank_;}

  /** @brief Clock error.
  * error = clock time - system time [s] */
  Double clockError(UInt idEpoch) const {return clk.at(idEpoch);}

  /** @brief set clock error.
  * error = observed clock time - system time [s] */
  void updateClockError(UInt idEpoch, Double deltaClock) {clk.at(idEpoch) += deltaClock;}

  /** @brief Estimated system time of received signals.
  * (receiver clock time - clock error) */
  Time timeCorrected(UInt idEpoch) const {return times.at(idEpoch) - seconds2time(clockError(idEpoch));}

  /** @brief Positions, velocities and orientations are given in the terrestrial reference frame. */
  Bool isEarthFixed() const {return isEarthFixed_;}

  /** @brief antenna reference point in TRF or CRF. */
  Vector3d position(UInt idEpoch) const {return pos.at(idEpoch) + global2local.at(idEpoch).inverseTransform(offset.at(idEpoch));}

  /** @brief velocity in TRF or CRF [m/s]. */
  Vector3d velocity(UInt idEpoch) const {return vel.at(idEpoch);}

  /** @brief Rotation from terrestrial reference frame (TRF) or celestial reference frame (CRF) to local horizont system (north, east, up). */
  Transform3d global2localFrame(UInt idEpoch) const {return global2local.at(idEpoch);}

  /** @brief Rotation from terrestrial reference frame (TRF) or celestial reference frame (CRF) to left-handed antenna system (usually north, east, up). */
  Transform3d global2antennaFrame(UInt idEpoch) const {return global2antenna.at(idEpoch);}

  /** @brief Transformation matrix for observed (composed) types from original transmitted types.
  * E.g. C2DG = C1CG - C1WG + C2WG.
  * Returns the @a typesTrans and the transformation matrix @a A (dimension: types.size() times typesTrans.size()). */
  virtual void signalComposition(UInt /*idEpoch*/, const std::vector<GnssType> &types, std::vector<GnssType> &typesTrans, Matrix &A) const;

  /** @brief All observations between receiver and transmitter at one epoch. */
  GnssObservation *observation(UInt idTrans, UInt idEpoch) const;

  /** @brief Max. observed epoch id +1. */
  UInt idEpochSize() const {return observations_.size();}

  /** @brief Max. observed transmitter id+1 at @a idEpoch. */
  UInt idTransmitterSize(UInt idEpoch) const {return (idEpoch < observations_.size()) ? observations_.at(idEpoch).size() : 0;}

  /** @brief Delete observation. */
  void deleteObservation(UInt idTrans, UInt idEpoch);

  /** @brief Delete all empty tracks (and ambiguities). */
  void deleteEmptyTracks();

  // Preprocessing
  // -------------
  /** @brief Compute observation equations (reduced observations). */
  class ObservationEquationList
  {
    std::vector<std::vector<std::unique_ptr<GnssObservationEquation>>> eqn;

  public:
    ObservationEquationList() {}
    ObservationEquationList(const GnssReceiver &receiver, const std::vector<GnssTransmitterPtr> &transmitters,
                            const std::function<Rotary3d(const Time &time)> &rotationCrf2Trf,
                            const std::function<void(GnssObservationEquation &eqn)> &reduceModels, GnssObservation::Group group);

    GnssObservationEquation *operator()(UInt idTrans, UInt idEpoch) const;
  };

  void preprocessingInfo(const std::string &info, UInt countEpochs=NULLINDEX, UInt countObservations=NULLINDEX, UInt countTracks=NULLINDEX);



  /** @brief Reads observations from a file. Member variable @a times must be set.
  * Initializes observations. Receiver and Transmitter positions, orientations, ... must be initialized beforehand.
  * Delete observations that don't match the types from receiver and transmitter definition. */
  void readObservations(const FileName &fileName, const std::vector<GnssTransmitterPtr> &transmitters, const std::function<Rotary3d(const Time &time)> &rotationCrf2Trf,
                        const Time &timeMargin, Angle elevationCutOff, const std::vector<GnssType> &useType, const std::vector<GnssType> &ignoreType, GnssObservation::Group group);

  /** @brief Substitute tracking attribute of signal type depending on receiver type */
  GnssType substituteSignal(GnssType type) const;

  /** @brief Simulate observations. Member variable @a times must be set.
  * Receiver and Transmitter positions, orientations, ... must be initialized beforehand.
  * Delete observations that don't match the types from receiver and transmitter definition. */
  void simulateObservations(const std::vector<GnssType> &types, Bool substituteTrackingMode,
                            NoiseGeneratorPtr noiseClock, NoiseGeneratorPtr noiseObs,
                            const std::vector<GnssTransmitterPtr> &transmitters,
                            const std::function<Rotary3d(const Time &time)> &rotationCrf2Trf,
                            const std::function<void(GnssObservationEquation &eqn)> &reduceModels,
                            UInt minObsCountPerTrack, Angle elevationCutOff, Angle elevationTrackMinimum,
                            const std::vector<GnssType> &useType,
                            const std::vector<GnssType> &ignoreType,
                            GnssObservation::Group group);

  /** @brief Estimate coarse receiver clock errors from a Precise Point Positioning (PPP) code solution.
  * If @p estimateKinematicPosition is TRUE, the receiver position is estimated at each epoch, otherwise it is estimated once for all epochs.*/
  std::vector<Vector3d> estimateInitialClockErrorFromCodeObservations(const std::vector<GnssTransmitterPtr> &transmitters,
                                                                      const std::function<Rotary3d(const Time &time)> &rotationCrf2Trf,
                                                                      const std::function<void(GnssObservationEquation &eqn)> &reduceModels,
                                                                      Double huber, Double huberPower, Bool estimateKinematicPosition);

  /** @brief Disable epochs if reduced code observations @p eqn exceed @p threshold (e.g. 100+ km for code cycle slips).
  * Disable receiver at epoch if @p outlierRatio or more of the observed satellites have gross code outliers. */
  void disableEpochsWithGrossCodeObservationOutliers(ObservationEquationList &eqn, Double threshold, Double outlierRatio=0.5);

  /** @brief Create tracks with continuously identical phase observations.
  * Tracks may contain gaps but must contain observations in at least @p minObsCountPerTrack epochs.
  * Extra types are included (e.g. L5*G), but tracks must have at least two others phases at different frequencies.*/
  void createTracks(const std::vector<GnssTransmitterPtr> &transmitters, UInt minObsCountPerTrack, const std::vector<GnssType> &extraTypes={});

  /** @brief delete track and all related observations. */
  void deleteTrack(UInt idTrack);

  /** @brief Removes tracks that never exceed @p minElevation (in radian). */
  void removeLowElevationTracks(ObservationEquationList &eqn, Angle minElevation);

  /** @brief Split a @p track at @p idEpochSplit into two new tracks.
  * Shortens the old track and returns the new track. Updates observation and observation equation @p eqn track assignments. */
  GnssTrackPtr splitTrack(ObservationEquationList &eqn, GnssTrackPtr track, UInt idEpochSplit);

  void linearCombinations(ObservationEquationList &eqnList, GnssTrackPtr track, const std::vector<GnssType> &extraTypes,
                          std::vector<GnssType> &typesPhase, std::vector<UInt> &idEpochs, Matrix &combinations, Double &cycles2tecu) const;

  void rangeAndTec(ObservationEquationList &eqnList, UInt idTrans, const std::vector<UInt> &idEpochs,
                   const std::vector<GnssType> &typesPhase, Vector &range, Vector &tec) const;

  void writeTracks(const FileName &fileName, ObservationEquationList &eqnList, const std::vector<GnssType> &extraTypes) const;

  /** @brief Splits tracks at detected cycle slips.
  * Based on all Melbourne-Wuebbena like combinations.
  * Solves the total variation regularized least-squares problem (total variation denoising).
  * In a second step smoothness of TEC is evaluated using a moving window peak/outlier detection based on AR model residuals.
  * @param eqnList Observation equations (reduced observations).
  * @param minObsCountPerTrack Removes tracks that don't have enough observations.
  * @param lambda Regularization parameter (@see @a totalVariationDenoising) (e.g. @p lambda = 5 for GPS ground stations).
  * @param windowSize Size of the moving window used for the TEC smoothness evaluation. If 0, TEC is not analyzed.
  * @param tecSigmaFactor Factor applied to moving standard deviation of AR model residuals to determine threshold for peak/outlier detection.
  * @param extraTypes GPS L5 observations are handled separately due to temporal changing bias.*/
  void cycleSlipsDetection(ObservationEquationList &eqnList, UInt minObsCountPerTrack, Double lambda, UInt windowSize, Double tecSigmaFactor, const std::vector<GnssType> &extraTypes={});
  void cycleSlipsDetection(ObservationEquationList &eqnList, GnssTrackPtr track, Double lambda, UInt windowSize, Double tecSigmaFactor, const std::vector<GnssType> &extraTypes);

  /** @brief repair cycle slip differences at same frequencies (e.g. between L1CG and L1WG).
  * Allows to reduce the number of integer ambiguities.
  * @param eqnList Observation equations (reduced observations). */
  void cycleSlipsRepairAtSameFrequency(ObservationEquationList &eqnList);

  /** @brief Track outlier detection based on robust least squares estimation.
  * Downweights outliers and prereduces observations by estimated integer ambiguities. */
  void trackOutlierDetection(const ObservationEquationList &eqn, const std::vector<GnssType> &ignoreTypes, Double huber, Double huberPower);

  /** @brief Robust least squares estimation for multiple epochs in a single system of equations. */
  static Matrix robustLeastSquares(const std::vector<Matrix> &A, const std::vector<Matrix> &l, const std::vector<UInt> &obsCount,
                                   Double huber, Double huberPower, UInt maxIter, Vector &sigma);

  /** @brief Total variation denoising.
  * Solves the total variation regularized least-squares problem.
  * Laurent Condat. A Direct Algorithm for 1D Total Variation Denoising. IEEE Signal Processing Letters,
  * Institute of Electrical and Electronics Engineers, 2013, 20 (11), pp.1054-1057. DOI: 10.1109/LSP.2013.2278339.
  * @param y noisy input time series.
  * @param lambda Regularization parameter.
  * @return denoised time series. */
  static Matrix totalVariationDenoising(const_MatrixSliceRef y, Double lambda);
}; //class GnssReceiver

/***** CLASS ***********************************/

class GnssAmbiguity;

/** @brief Class for tracking phases along epochs. */
class GnssTrack
{
public:
  GnssReceiver         *receiver;
  GnssTransmitter      *transmitter;
  UInt                  idEpochStart, idEpochEnd;
  std::vector<GnssType> types;
  GnssAmbiguity        *ambiguity;

  GnssTrack(GnssReceiver *_receiver, GnssTransmitter *_transmitter, UInt _idEpochStart, UInt _idEpochEnd, const std::vector<GnssType> &_types);
 ~GnssTrack();

  UInt countObservations() const;
  void removeAmbiguitiesFromObservations(const std::vector<GnssType> &type, const std::vector<Double> &value);
};

/***** CLASS ***********************************/

/** @brief Abstract class for ambiguities.
* @a track is owns this class and delete ambiguity when destructed. */
class GnssAmbiguity
{
public:
  GnssTrack *track;

  explicit GnssAmbiguity(GnssTrack *track) : track(track) {track->ambiguity = this;}
  virtual ~GnssAmbiguity() {}
  virtual Vector ambiguities(const std::vector<GnssType> &type) const = 0;
};

/***********************************************/

/// @}

#endif /* __GROOPS___ */
