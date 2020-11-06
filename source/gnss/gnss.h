/***********************************************/
/**
* @file gnss.h
*
* @brief global navigation satellite system classes.
*
* @author Torsten Mayer-Guerr
* @date 2010-08-03
*
*/
/***********************************************/

#ifndef __GROOPS_GNSS__
#define __GROOPS_GNSS__

#include "parallel/parallel.h"
#include "base/gnssType.h"
#include "base/parameterName.h"

/** @addtogroup gnssGroup */
/// @{

/***** TYPES ***********************************/

class Gnss;
typedef std::shared_ptr<Gnss> GnssPtr;

class GnssParametrizationIonosphere;
typedef std::shared_ptr<GnssParametrizationIonosphere> GnssParametrizationIonospherePtr;

class GnssParametrizationEarthRotation;
typedef std::shared_ptr<GnssParametrizationEarthRotation> GnssParametrizationEarthRotationPtr;

class GnssParametrizationGravityField;
typedef std::shared_ptr<GnssParametrizationGravityField> GnssParametrizationGravityFieldPtr;

class GnssParametrizationAmbiguities;
typedef std::shared_ptr<GnssParametrizationAmbiguities> GnssParametrizationAmbiguitiesPtr;

class MatrixDistributed;

/***** CLASS ***********************************/

class Gnss
{
public:
  class SingleObservation;
  class Observation;
  class ObservationEquation;
  class Receiver;
  class Transmitter;
  class Track;
  class Ambiguity;
  class DesignMatrix;
  class NormalEquationInfo;
  class Parametrization;

  typedef std::shared_ptr<Receiver>        ReceiverPtr;
  typedef std::shared_ptr<Transmitter>     TransmitterPtr;
  typedef std::shared_ptr<Parametrization> ParametrizationPtr;

  typedef UInt AnalysisType;
  constexpr static AnalysisType ANALYSIS_CODE  = 1<<0; // use code observations only
  constexpr static AnalysisType ANALYSIS_PHASE = 1<<1;
  constexpr static AnalysisType ANALYSIS_IONO  = 1<<2;
  constexpr static AnalysisType ANALYSIS_RAW   = ~0;   // all bits set: use all available observations as they are

  // ======================

  class ParameterIndex
  {
    UInt index;
  public:
    ParameterIndex(UInt idx=NULLINDEX) : index(idx) {}
    explicit operator bool() const {return index != NULLINDEX;}

    friend class NormalEquationInfo;
    friend class DesignMatrix;
  };

  // ======================

  class NormalEquationInfo
  {
  public:
    typedef UInt64 EstimationType;
    static std::vector<std::tuple<Gnss::NormalEquationInfo::EstimationType, std::string, std::string, std::string>> estimationTypes;

    // ATTENTION: if you add a new EstimationType here, also an entry in Gnss::NormalEquationInfo::estimationTypes (in gnss.cpp)
    constexpr static EstimationType ESTIMATE_RECEIVER_CLOCK                         = EstimationType(1)<<0;
    constexpr static EstimationType ESTIMATE_RECEIVER_KINEMATICPOSITION             = EstimationType(1)<<1;
    constexpr static EstimationType ESTIMATE_RECEIVER_POSITION                      = EstimationType(1)<<2;
    constexpr static EstimationType ESTIMATE_RECEIVER_TROPOSPHERE_WET               = EstimationType(1)<<3;
    constexpr static EstimationType ESTIMATE_RECEIVER_TROPOSPHERE_GRADIENT          = EstimationType(1)<<4;
    constexpr static EstimationType ESTIMATE_RECEIVER_SIGNALBIAS                    = EstimationType(1)<<5;
    constexpr static EstimationType ESTIMATE_RECEIVER_TECBIAS                       = EstimationType(1)<<6;
    constexpr static EstimationType ESTIMATE_RECEIVER_ANTENNACENTERVARIATIONS       = EstimationType(1)<<7;
    constexpr static EstimationType ESTIMATE_TRANSMITTER_CLOCK                      = EstimationType(1)<<8;
    constexpr static EstimationType ESTIMATE_TRANSMITTER_KINEMATICPOSITION          = EstimationType(1)<<9;
    constexpr static EstimationType ESTIMATE_TRANSMITTER_ORBIT                      = EstimationType(1)<<10;
    constexpr static EstimationType ESTIMATE_TRANSMITTER_CLOCKMODEL                 = EstimationType(1)<<11;
    constexpr static EstimationType ESTIMATE_TRANSMITTER_SIGNALBIAS                 = EstimationType(1)<<12;
    constexpr static EstimationType ESTIMATE_TRANSMITTER_TECBIAS                    = EstimationType(1)<<13;
    constexpr static EstimationType ESTIMATE_TRANSMITTER_SIGNALBIASMODEL            = EstimationType(1)<<14;
    constexpr static EstimationType ESTIMATE_TRANSMITTER_ANTENNACENTERVARIATIONS    = EstimationType(1)<<15;
    constexpr static EstimationType ESTIMATE_EARTHROTATION_POLE                     = EstimationType(1)<<16;
    constexpr static EstimationType ESTIMATE_EARTHROTATION_UT1                      = EstimationType(1)<<17;
    constexpr static EstimationType ESTIMATE_EARTHROTATION_NUTATION                 = EstimationType(1)<<18;
    constexpr static EstimationType ESTIMATE_GRAVITY                                = EstimationType(1)<<19;
    constexpr static EstimationType ESTIMATE_AMBIGUITIES                            = EstimationType(1)<<20;
    constexpr static EstimationType ESTIMATE_IONOSPHERE_STEC                        = EstimationType(1)<<21;
    constexpr static EstimationType ESTIMATE_IONOSPHERE_VTEC                        = EstimationType(1)<<22;
    constexpr static EstimationType ESTIMATE_IONOSPHERE_TECMAP                      = EstimationType(1)<<23;
    constexpr static EstimationType CONSTRAINT_RECEIVER_SIGNALBIAS_DIFFCLOCK        = EstimationType(1)<<24;
    constexpr static EstimationType CONSTRAINT_RECEIVER_POSITION_NONETTRANSLATION   = EstimationType(1)<<25;
    constexpr static EstimationType CONSTRAINT_RECEIVER_POSITION_NONETROTATION      = EstimationType(1)<<26;
    constexpr static EstimationType CONSTRAINT_TRANSMITTER_CLOCK_ZEROMEAN           = EstimationType(1)<<27;
    constexpr static EstimationType CONSTRAINT_TRANSMITTER_CLOCKMODEL               = EstimationType(1)<<28;
    constexpr static EstimationType CONSTRAINT_TRANSMITTER_SIGNALBIAS_ZEROMEAN      = EstimationType(1)<<29;
    constexpr static EstimationType CONSTRAINT_TRANSMITTER_SIGNALBIASMODEL_ZEROMEAN = EstimationType(1)<<30;
    constexpr static EstimationType CONSTRAINT_IONOSPHERE_STEC                      = EstimationType(1)<<31;
    constexpr static EstimationType CONSTRAINT_OTHER                                = EstimationType(1)<<32;

    constexpr static EstimationType MASK_ALL            = ~0;
    constexpr static EstimationType MASK_RECEIVER       = ESTIMATE_RECEIVER_CLOCK | ESTIMATE_RECEIVER_KINEMATICPOSITION | ESTIMATE_RECEIVER_POSITION
                                                        | ESTIMATE_RECEIVER_TROPOSPHERE_WET | ESTIMATE_RECEIVER_TROPOSPHERE_GRADIENT
                                                        | ESTIMATE_RECEIVER_SIGNALBIAS | ESTIMATE_RECEIVER_TECBIAS
                                                        | ESTIMATE_RECEIVER_ANTENNACENTERVARIATIONS
                                                        | CONSTRAINT_RECEIVER_SIGNALBIAS_DIFFCLOCK
                                                        | CONSTRAINT_RECEIVER_POSITION_NONETTRANSLATION | CONSTRAINT_RECEIVER_POSITION_NONETROTATION;
    constexpr static EstimationType MASK_TRANSMITTER    = ESTIMATE_TRANSMITTER_CLOCK | ESTIMATE_TRANSMITTER_KINEMATICPOSITION | ESTIMATE_TRANSMITTER_ORBIT| ESTIMATE_TRANSMITTER_CLOCKMODEL
                                                        | ESTIMATE_TRANSMITTER_SIGNALBIAS | ESTIMATE_TRANSMITTER_TECBIAS | ESTIMATE_TRANSMITTER_SIGNALBIASMODEL
                                                        | ESTIMATE_TRANSMITTER_ANTENNACENTERVARIATIONS
                                                        | CONSTRAINT_TRANSMITTER_CLOCK_ZEROMEAN | CONSTRAINT_TRANSMITTER_CLOCKMODEL
                                                        | CONSTRAINT_TRANSMITTER_SIGNALBIAS_ZEROMEAN | CONSTRAINT_TRANSMITTER_SIGNALBIASMODEL_ZEROMEAN;
    constexpr static EstimationType MASK_EPOCH          = ESTIMATE_RECEIVER_CLOCK    | ESTIMATE_RECEIVER_KINEMATICPOSITION
                                                        | ESTIMATE_TRANSMITTER_CLOCK | ESTIMATE_TRANSMITTER_KINEMATICPOSITION
                                                        | CONSTRAINT_TRANSMITTER_CLOCK_ZEROMEAN | CONSTRAINT_TRANSMITTER_CLOCKMODEL;
    constexpr static EstimationType MASK_CONSTRAINT     = CONSTRAINT_RECEIVER_SIGNALBIAS_DIFFCLOCK
                                                        | CONSTRAINT_RECEIVER_POSITION_NONETTRANSLATION | CONSTRAINT_RECEIVER_POSITION_NONETROTATION
                                                        | CONSTRAINT_TRANSMITTER_CLOCK_ZEROMEAN | CONSTRAINT_TRANSMITTER_CLOCKMODEL
                                                        | CONSTRAINT_TRANSMITTER_SIGNALBIAS_ZEROMEAN | CONSTRAINT_TRANSMITTER_SIGNALBIASMODEL_ZEROMEAN
                                                        | CONSTRAINT_IONOSPHERE_STEC | CONSTRAINT_OTHER;
    constexpr static EstimationType MASK_SINGLERECEIVER = (MASK_RECEIVER | ESTIMATE_AMBIGUITIES | ESTIMATE_IONOSPHERE_STEC | ESTIMATE_IONOSPHERE_VTEC | MASK_CONSTRAINT)
                                                        & ~CONSTRAINT_RECEIVER_SIGNALBIAS_DIFFCLOCK & ~CONSTRAINT_RECEIVER_POSITION_NONETTRANSLATION & ~CONSTRAINT_RECEIVER_POSITION_NONETROTATION;

    Parallel::CommunicatorPtr comm;
    AnalysisType              analysisType;
    EstimationType            estimationType;
    std::vector<Byte>         estimateReceiver;  // subset of observations
    std::vector<UInt>         idEpochs;          // epochs to estimate
    UInt                      defaultBlockSizeEpoch;
    UInt                      defaultBlockSizeInterval;
    UInt                      defaultBlockSizeAmbiguity;
    UInt                      defaultBlockReceiverCount;
    UInt                      defaultBlockCountReduction;
    Bool                      keepEpochNormalsInMemory;
    Bool                      accumulateEpochObservations;

    NormalEquationInfo(UInt countEpoch, UInt countReceiver, UInt countTransmitter, Parallel::CommunicatorPtr comm);

    void initNewParameterNames();
    ParameterIndex parameterNamesEpochReceiver   (UInt idEpoch, UInt idRecv,  const std::vector<ParameterName> &parameterNames)              {return addParameters(idEpoch,   idRecv,    NULLINDEX, parameterNames);}
    ParameterIndex parameterNamesEpochTransmitter(UInt idEpoch, UInt idTrans, const std::vector<ParameterName> &parameterNames)              {return addParameters(idEpoch,   NULLINDEX, idTrans,   parameterNames);}
    ParameterIndex parameterNamesEpochOther      (UInt idEpoch, const std::vector<ParameterName> &parameterNames)                            {return addParameters(idEpoch,   NULLINDEX, NULLINDEX, parameterNames);}
    ParameterIndex parameterNamesReceiver        (UInt idRecv,  const std::vector<ParameterName> &parameterNames)                            {return addParameters(NULLINDEX, idRecv,    NULLINDEX, parameterNames);}
    ParameterIndex parameterNamesTransmitter     (UInt idTrans, const std::vector<ParameterName> &parameterNames)                            {return addParameters(NULLINDEX, NULLINDEX, idTrans,   parameterNames);}
    ParameterIndex parameterNamesOther           (const std::vector<ParameterName> &parameterNames)                                          {return addParameters(NULLINDEX, NULLINDEX, NULLINDEX, parameterNames);}
    ParameterIndex parameterNamesAmbiguity       (UInt idEpoch, UInt idRecv, UInt idTrans, const std::vector<ParameterName> &parameterNames) {return addParameters(idEpoch,   idRecv,    idTrans,   parameterNames);}
    void calculateIndex();

    UInt block(const ParameterIndex &index) const {return block_.at(index.index);}
    UInt index(const ParameterIndex &index) const {return index_.at(index.index);}
    UInt count(const ParameterIndex &index) const {return count_.at(index.index);}

    const std::vector<ParameterName> &parameterNames() const {return parameterNames_;}
    const std::vector<UInt>          &blockIndices()   const {return blockIndices_;}

    UInt parameterCount()   const {return blockIndices_.back();}                      //!< Number of rows/columns (dimension) of distributed matrix
    UInt blockIndex(UInt i) const {return blockIndices_.at(i);}                       //!< Start index of block @a i
    UInt blockSize(UInt i)  const {return blockIndices_.at(i+1)-blockIndices_.at(i);} //!< Size of block @a i
    UInt blockCount()       const {return blockIndices_.size()-1;}                    //!< Number of block rows/columns

    UInt blockCountEpoch(UInt idEpoch) const {return blockCountEpoch_.at(idEpoch);}
    UInt blockInterval()  const {return blockInterval_;}
    UInt blockAmbiguity() const {return blockAmbiguity_;}

    void removeEpochParameters();

  private:
    std::vector<std::tuple<UInt, UInt, UInt, UInt, std::vector<ParameterName>>> parameters; // idEpoch, idRecv, idTrans, idx, name
    std::vector<UInt>          block_, index_, count_;
    std::vector<ParameterName> parameterNames_;
    std::vector<UInt>          blockIndices_;
    std::vector<UInt>          blockCountEpoch_;
    UInt                       blockInterval_, blockAmbiguity_;
    UInt                       countTransmitter_;

    ParameterIndex addParameters(UInt idEpoch, UInt idRecv, UInt idTrans, const std::vector<ParameterName> &parameterNames);
  };

  // ======================

  class Parametrization
  {
    friend class Gnss;
    Gnss   *_gnss; // is set by Gnss::registerParametrization

    public:
    virtual ~Parametrization() {}

    /** @brief Reference to the complete GNSS system.
    * After registerParametrization. */
    Gnss &gnss() const;

    // is called by Gnss::initInterval by every node in comm
    // several steps: first early is called, transmitter are initialized secondly, then receivers (with observations) and last all others
    virtual void initIntervalEarly      (AnalysisType analysisType, const std::vector<Time> &times, const Time &timeMargin, Parallel::CommunicatorPtr comm);
    virtual void initIntervalTransmitter(AnalysisType analysisType, const std::vector<Time> &times, const Time &timeMargin, Parallel::CommunicatorPtr comm);
    virtual void initIntervalReceiver   (AnalysisType analysisType, const std::vector<Time> &times, const Time &timeMargin, Parallel::CommunicatorPtr comm);
    virtual Bool initIntervalDisabling  (AnalysisType analysisType, const std::vector<Time> &times, const Time &timeMargin, Parallel::CommunicatorPtr comm);
    virtual void initIntervalLate       (AnalysisType analysisType, const std::vector<Time> &times, const Time &timeMargin, Parallel::CommunicatorPtr comm);

    // called by every node in comm of normalEquationInfo
    // comm might be a subset of comm from initInterval (for single receiver)
    // returned normalEquationInfo must be the same on every node
    virtual void initParameter(NormalEquationInfo &normalEquationInfo);

    // called by every node in comm of normalEquationInfo
    // @param normalEquationInfo from @a initParameter
    // @param[out] x0 distributed at nodes, @a Parallel::reduceSum must be called afterwards
    virtual void aprioriParameter(const NormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const;

    // called by only one node in comm of normalEquationInfo
    virtual Bool isDesignMatrix(const NormalEquationInfo &normalEquationInfo, UInt idRecv, UInt idTrans, UInt idEpoch) const;

    // called by only one node in comm of normalEquationInfo
    virtual void designMatrix(const NormalEquationInfo &normalEquationInfo, const ObservationEquation &eqn, DesignMatrix &A) const;

    // called by every node in comm of normalEquationInfo
    // returned each equation only at one arbritary node
    virtual void observationEquationEpoch(const NormalEquationInfo &normalEquationInfo, UInt idEpoch, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const;

    // called by every node in comm of normalEquationInfo
    // returned each equation only at one arbritary node
    virtual void observationEquation(const NormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const;

    // called by every node in comm of normalEquationInfo
    // @return max absolute parameter change in meter (for convergence tracking)
    virtual Double updateParameter(const NormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz, Bool printStatistics);

    // called by every node in comm of normalEquationInfo
    virtual void updateCovariance(const NormalEquationInfo &normalEquationInfo, const MatrixDistributed &covariance);

    // save estimated parameters
    // called by every node in comm of normalEquationInfo
    virtual void writeResults(const NormalEquationInfo &normalEquationInfo, const std::string &suffix);
  };

  // ===========================================================

  std::vector<Time>                   times;           // Epochs
  std::vector<ReceiverPtr>            receiver;        // stations & LEOs
  std::vector<TransmitterPtr>         transmitter;     // GNSS satellites
  std::vector<ParametrizationPtr>     parametrization; // parameters
  std::vector<std::vector<std::vector<GnssType>>> typesRecvTrans; // for each receiver and transmitter: used types (receiver types)
  GnssParametrizationIonospherePtr    ionosphere;
  GnssParametrizationEarthRotationPtr earthRotation;
  GnssParametrizationGravityFieldPtr  gravityField;
  GnssParametrizationAmbiguitiesPtr   ambiguities;     // ambiguity parameters

  Gnss() {}
 ~Gnss() {}

  // Parametrization
  // ---------------
  void   init(const std::vector<ReceiverPtr> &receiver, const std::vector<TransmitterPtr> &transmitter,
              const std::vector<ParametrizationPtr> &parametrization,
              GnssParametrizationIonospherePtr ionosphere,
              GnssParametrizationEarthRotationPtr gnssEarthRotation,
              GnssParametrizationGravityFieldPtr gnssGravityField,
              GnssParametrizationAmbiguitiesPtr gnssAmbiguities);
  void   initInterval             (AnalysisType analysisType, const std::vector<Time> &times, const Time &timeMargin, Parallel::CommunicatorPtr comm);
  void   initParameter            (NormalEquationInfo &normalEquationInfo);
  // called by every node in comm of normalEquationInfo
  // @return only valid at master
  Vector aprioriParameter         (const NormalEquationInfo &normalEquationInfo) const;
  Bool   isDesignMatrix           (const NormalEquationInfo &normalEquationInfo, UInt idRecv, UInt idTrans, UInt idEpoch) const;
  Bool   basicObservationEquations(const NormalEquationInfo &normalEquationInfo, UInt idRecv, UInt idTrans, UInt idEpoch, Gnss::ObservationEquation &eqn) const;
  void   designMatrix             (const NormalEquationInfo &normalEquationInfo, const ObservationEquation &eqn, DesignMatrix &A) const;
  void   observationEquationEpoch (const NormalEquationInfo &normalEquationInfo, UInt idEpoch, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const;
  void   observationEquation      (const NormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const;
  Double ambiguityResolve         (const NormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount, Bool dryRun=FALSE);
  Double updateParameter          (const NormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz, Bool printStatistics);
  void   updateCovariance         (const NormalEquationInfo &normalEquationInfo, const MatrixDistributed &covariance);
  void   writeResults             (const NormalEquationInfo &normalEquationInfo, const std::string &suffix="");

  // ===========================================================

  /** @brief sorted list of used types. */
  std::vector<GnssType> types(const GnssType mask=GnssType::ALL) const;

  /** @brief Replaces observed (composed) types by orignal transmitted types.
  * E.g. C2DG = C1CG - C1WG + C2CW.
  * Phase types returned without tracking attribute. */
  static std::vector<GnssType> replaceCompositeSignals(const std::vector<GnssType> &types);

  /** @brief Transformation matrix for observed (composed) types from orignal transmitted types.
  * E.g. C2DG = C1CG - C1WG + C2WG.
  * Returns the @a typesTrans and the transformation matrix @a A (dimension: types.size() times typesTrans.size()).
  * Assumes equally power of the two orginal signals for the composed C*X* signals. */
  static void defaultSignalComposition(const std::vector<GnssType> &types, std::vector<GnssType> &typesTrans, Matrix &A);

  /** @brief Groups the observed/transmitted types into independent systems without common types.
  * First example: A receiver @a idDevice=idRecv() (@a deviceIsReceiver=TRUE) observes GPS, GALIEO, GLONASS satellites.
  * Second example: A Gallieo transmitter @a idDevice=idTrans() (@a deviceIsReceiver=FALSE) is observed by a group of stations with L1 and L5 and by another group with L1 and L8 signals.
  * @a typesRecvTrans are the obsereved types for each receiver and transmitter. */
  std::vector<std::vector<GnssType>> groupTypes(const std::vector<std::vector<std::vector<GnssType>>> &typesRecvTrans, UInt idDevice, Bool deviceIsReceiver) const;

  void signalBiasParameter(const std::vector<GnssType> &obsTypes, const std::vector<std::vector<GnssType>> &groupTypes,
                           Bool eliminateClock, Bool estimateCodeBias, Bool estimatePhaseBias, Bool estimateTecBias,
                           Matrix &A, std::vector<ParameterName> &parameterName) const;

  static void checkMaxChange(Double &maxChange, std::string &info, Bool printStatistics, Parallel::CommunicatorPtr comm);
};

/// @}

/***********************************************/
/*** INLINES ***********************************/
/***********************************************/

inline void   Gnss::Parametrization::initIntervalEarly       (Gnss::AnalysisType /*analysisType*/, const std::vector<Time> &/*times*/, const Time &/*timeMargin*/, Parallel::CommunicatorPtr /*comm*/) {}
inline void   Gnss::Parametrization::initIntervalTransmitter (Gnss::AnalysisType /*analysisType*/, const std::vector<Time> &/*times*/, const Time &/*timeMargin*/, Parallel::CommunicatorPtr /*comm*/) {}
inline void   Gnss::Parametrization::initIntervalReceiver    (Gnss::AnalysisType /*analysisType*/, const std::vector<Time> &/*times*/, const Time &/*timeMargin*/, Parallel::CommunicatorPtr /*comm*/) {}
inline Bool   Gnss::Parametrization::initIntervalDisabling   (Gnss::AnalysisType /*analysisType*/, const std::vector<Time> &/*times*/, const Time &/*timeMargin*/, Parallel::CommunicatorPtr /*comm*/) {return FALSE;}
inline void   Gnss::Parametrization::initIntervalLate        (Gnss::AnalysisType /*analysisType*/, const std::vector<Time> &/*times*/, const Time &/*timeMargin*/, Parallel::CommunicatorPtr /*comm*/) {}
inline void   Gnss::Parametrization::initParameter           (NormalEquationInfo &/*normalEquationInfo*/) {}
inline void   Gnss::Parametrization::aprioriParameter        (const NormalEquationInfo &/*normalEquationInfo*/, MatrixSliceRef /*x0*/) const {}
inline Bool   Gnss::Parametrization::isDesignMatrix          (const NormalEquationInfo &/*normalEquationInfo*/, UInt /*idRecv*/, UInt /*idTrans*/, UInt /*idEpoch*/) const {return FALSE;}
inline void   Gnss::Parametrization::designMatrix            (const NormalEquationInfo &/*normalEquationInfo*/, const ObservationEquation &/*eqn*/, DesignMatrix &/*A*/) const {}
inline void   Gnss::Parametrization::observationEquationEpoch(const NormalEquationInfo &/*normalEquationInfo*/, UInt /*idEpoch*/, MatrixDistributed &/*normals*/, std::vector<Matrix> &/*n*/, Double &/*lPl*/, UInt &/*obsCount*/) const {}
inline void   Gnss::Parametrization::observationEquation     (const NormalEquationInfo &/*normalEquationInfo*/, MatrixDistributed &/*normals*/, std::vector<Matrix> &/*n*/, Double &/*lPl*/, UInt &/*obsCount*/) const {}
inline Double Gnss::Parametrization::updateParameter         (const NormalEquationInfo &/*normalEquationInfo*/, const_MatrixSliceRef /*x*/, const_MatrixSliceRef /*Wz*/, Bool /*printStatistics*/) {return 0;}
inline void   Gnss::Parametrization::updateCovariance        (const NormalEquationInfo &/*normalEquationInfo*/, const MatrixDistributed &/*covariance*/) {}
inline void   Gnss::Parametrization::writeResults            (const NormalEquationInfo &/*normalEquationInfo*/, const std::string &/*suffix*/) {}

/***********************************************/

#endif /* __GROOPS___ */
