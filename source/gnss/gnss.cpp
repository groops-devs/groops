/***********************************************/
/**
* @file gnss.cpp
*
* @brief global navigation satellite system.
*
* @author Torsten Mayer-Guerr
* @date 2010-08-03
*
*/
/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "inputOutput/logging.h"
#include "gnss.h"
#include "gnss/gnssObservation.h"
#include "gnss/gnssTransmitter.h"
#include "gnss/gnssReceiver.h"
#include "gnss/gnssParametrizationIonosphere.h"
#include "gnss/gnssParametrizationEarthRotation.h"
#include "gnss/gnssParametrizationGravityField.h"
#include "gnss/gnssParametrizationAmbiguities.h"

/***********************************************/

constexpr Gnss::NormalEquationInfo::EstimationType Gnss::NormalEquationInfo::ESTIMATE_RECEIVER_CLOCK;
constexpr Gnss::NormalEquationInfo::EstimationType Gnss::NormalEquationInfo::ESTIMATE_RECEIVER_KINEMATICPOSITION;
constexpr Gnss::NormalEquationInfo::EstimationType Gnss::NormalEquationInfo::ESTIMATE_RECEIVER_POSITION;
constexpr Gnss::NormalEquationInfo::EstimationType Gnss::NormalEquationInfo::CONSTRAINT_RECEIVER_POSITION_NONETTRANSLATION;
constexpr Gnss::NormalEquationInfo::EstimationType Gnss::NormalEquationInfo::CONSTRAINT_RECEIVER_POSITION_NONETROTATION;
constexpr Gnss::NormalEquationInfo::EstimationType Gnss::NormalEquationInfo::ESTIMATE_RECEIVER_TROPOSPHERE_WET;
constexpr Gnss::NormalEquationInfo::EstimationType Gnss::NormalEquationInfo::ESTIMATE_RECEIVER_TROPOSPHERE_GRADIENT;
constexpr Gnss::NormalEquationInfo::EstimationType Gnss::NormalEquationInfo::ESTIMATE_RECEIVER_SIGNALBIAS;
constexpr Gnss::NormalEquationInfo::EstimationType Gnss::NormalEquationInfo::CONSTRAINT_RECEIVER_SIGNALBIAS_DIFFCLOCK;
constexpr Gnss::NormalEquationInfo::EstimationType Gnss::NormalEquationInfo::ESTIMATE_RECEIVER_TECBIAS;
constexpr Gnss::NormalEquationInfo::EstimationType Gnss::NormalEquationInfo::ESTIMATE_RECEIVER_ANTENNACENTERVARIATIONS;
constexpr Gnss::NormalEquationInfo::EstimationType Gnss::NormalEquationInfo::ESTIMATE_TRANSMITTER_CLOCK;
constexpr Gnss::NormalEquationInfo::EstimationType Gnss::NormalEquationInfo::CONSTRAINT_TRANSMITTER_CLOCK_ZEROMEAN;
constexpr Gnss::NormalEquationInfo::EstimationType Gnss::NormalEquationInfo::CONSTRAINT_TRANSMITTER_CLOCKMODEL;
constexpr Gnss::NormalEquationInfo::EstimationType Gnss::NormalEquationInfo::ESTIMATE_TRANSMITTER_KINEMATICPOSITION;
constexpr Gnss::NormalEquationInfo::EstimationType Gnss::NormalEquationInfo::ESTIMATE_TRANSMITTER_ORBIT;
constexpr Gnss::NormalEquationInfo::EstimationType Gnss::NormalEquationInfo::ESTIMATE_TRANSMITTER_CLOCKMODEL;
constexpr Gnss::NormalEquationInfo::EstimationType Gnss::NormalEquationInfo::ESTIMATE_TRANSMITTER_SIGNALBIAS;
constexpr Gnss::NormalEquationInfo::EstimationType Gnss::NormalEquationInfo::CONSTRAINT_TRANSMITTER_SIGNALBIAS_ZEROMEAN;
constexpr Gnss::NormalEquationInfo::EstimationType Gnss::NormalEquationInfo::ESTIMATE_TRANSMITTER_TECBIAS;
constexpr Gnss::NormalEquationInfo::EstimationType Gnss::NormalEquationInfo::ESTIMATE_TRANSMITTER_SIGNALBIASMODEL;
constexpr Gnss::NormalEquationInfo::EstimationType Gnss::NormalEquationInfo::CONSTRAINT_TRANSMITTER_SIGNALBIASMODEL_ZEROMEAN;
constexpr Gnss::NormalEquationInfo::EstimationType Gnss::NormalEquationInfo::ESTIMATE_TRANSMITTER_ANTENNACENTERVARIATIONS;
constexpr Gnss::NormalEquationInfo::EstimationType Gnss::NormalEquationInfo::ESTIMATE_EARTHROTATION_POLE;
constexpr Gnss::NormalEquationInfo::EstimationType Gnss::NormalEquationInfo::ESTIMATE_EARTHROTATION_UT1;
constexpr Gnss::NormalEquationInfo::EstimationType Gnss::NormalEquationInfo::ESTIMATE_EARTHROTATION_NUTATION;
constexpr Gnss::NormalEquationInfo::EstimationType Gnss::NormalEquationInfo::ESTIMATE_GRAVITY;
constexpr Gnss::NormalEquationInfo::EstimationType Gnss::NormalEquationInfo::ESTIMATE_AMBIGUITIES;
constexpr Gnss::NormalEquationInfo::EstimationType Gnss::NormalEquationInfo::ESTIMATE_IONOSPHERE_STEC;
constexpr Gnss::NormalEquationInfo::EstimationType Gnss::NormalEquationInfo::CONSTRAINT_IONOSPHERE_STEC;
constexpr Gnss::NormalEquationInfo::EstimationType Gnss::NormalEquationInfo::ESTIMATE_IONOSPHERE_VTEC;
constexpr Gnss::NormalEquationInfo::EstimationType Gnss::NormalEquationInfo::ESTIMATE_IONOSPHERE_TECMAP;
constexpr Gnss::NormalEquationInfo::EstimationType Gnss::NormalEquationInfo::CONSTRAINT_OTHER;

/***********************************************/

std::vector<std::tuple<Gnss::NormalEquationInfo::EstimationType, std::string, std::string, std::string>> Gnss::NormalEquationInfo::estimationTypes = {
//                         flag,                                            configElement,         default, comment
// ----------------------------------------------------------------------------------------------------------------
{Gnss::NormalEquationInfo::ESTIMATE_RECEIVER_CLOCK,                         "estimateReceiverClock",                        "1", ""},
{Gnss::NormalEquationInfo::ESTIMATE_RECEIVER_KINEMATICPOSITION,             "estimateReceiverKinematicPosition",            "1", ""},
{Gnss::NormalEquationInfo::ESTIMATE_RECEIVER_POSITION,                      "estimateReceiverPosition",                     "1", ""},
{Gnss::NormalEquationInfo::CONSTRAINT_RECEIVER_POSITION_NONETTRANSLATION,   "constraintReceiverNoNetTranslation",           "1", ""},
{Gnss::NormalEquationInfo::CONSTRAINT_RECEIVER_POSITION_NONETROTATION,      "constraintReceiverNoNetRotation",              "1", ""},
{Gnss::NormalEquationInfo::ESTIMATE_RECEIVER_TROPOSPHERE_WET,               "estimateReceiverTroposphereWet",               "1", ""},
{Gnss::NormalEquationInfo::ESTIMATE_RECEIVER_TROPOSPHERE_GRADIENT,          "estimateReceiverTroposphereGradient",          "1", ""},
{Gnss::NormalEquationInfo::ESTIMATE_RECEIVER_SIGNALBIAS,                    "estimateReceiverSignalBias",                   "1", ""},
{Gnss::NormalEquationInfo::CONSTRAINT_RECEIVER_SIGNALBIAS_DIFFCLOCK,        "constraintReceiverSignalBiasDiffClock",        "1", "separate clock errors and biases of different satellite constellations"},
{Gnss::NormalEquationInfo::ESTIMATE_RECEIVER_TECBIAS,                       "estimateReceiverTecBias",                      "1", ""},
{Gnss::NormalEquationInfo::ESTIMATE_RECEIVER_ANTENNACENTERVARIATIONS,       "estimateReceiverAntennaCenterVariations",      "1", ""},
{Gnss::NormalEquationInfo::ESTIMATE_TRANSMITTER_CLOCK,                      "estimateTransmitterClock",                     "1", ""},
{Gnss::NormalEquationInfo::CONSTRAINT_TRANSMITTER_CLOCK_ZEROMEAN,           "constraintTransmitterClockZeroMean",           "1", "defines the time system"},
{Gnss::NormalEquationInfo::CONSTRAINT_TRANSMITTER_CLOCKMODEL,               "constraintTransmitterClockModel",              "1", ""},
{Gnss::NormalEquationInfo::ESTIMATE_TRANSMITTER_KINEMATICPOSITION,          "estimateTransmitterKinematicPosition",         "1", ""},
{Gnss::NormalEquationInfo::ESTIMATE_TRANSMITTER_ORBIT,                      "estimateTransmitterOrbit",                     "1", "incl. SRP parameters and pulses"},
{Gnss::NormalEquationInfo::ESTIMATE_TRANSMITTER_CLOCKMODEL,                 "estimateTransmitterClockModel",                "1", ""},
{Gnss::NormalEquationInfo::ESTIMATE_TRANSMITTER_SIGNALBIAS,                 "estimateTransmitterSignalBias",                "1", ""},
{Gnss::NormalEquationInfo::CONSTRAINT_TRANSMITTER_SIGNALBIAS_ZEROMEAN,      "constraintTransmitterSignalBiasZeroMean",      "1", "separate transmitter and receiver biases"},
{Gnss::NormalEquationInfo::ESTIMATE_TRANSMITTER_TECBIAS,                    "estimateTransmitterTecBias",                   "1", ""},
{Gnss::NormalEquationInfo::ESTIMATE_TRANSMITTER_SIGNALBIASMODEL,            "estimateTransmitterBiasModel",                 "1", ""},
{Gnss::NormalEquationInfo::CONSTRAINT_TRANSMITTER_SIGNALBIASMODEL_ZEROMEAN, "constraintTransmitterSignalBiasModelZeroMean", "1", "separate temporal biases from constant biases"},
{Gnss::NormalEquationInfo::ESTIMATE_TRANSMITTER_ANTENNACENTERVARIATIONS,    "estimateTransmitterAntennaCenterVariations",   "1", ""},
{Gnss::NormalEquationInfo::ESTIMATE_EARTHROTATION_POLE,                     "estimateEarthRotationPole",                    "1", ""},
{Gnss::NormalEquationInfo::ESTIMATE_EARTHROTATION_UT1,                      "estimateEarthRotationUT1",                     "1", ""},
{Gnss::NormalEquationInfo::ESTIMATE_EARTHROTATION_NUTATION,                 "estimateEarthRotationNutation",                "1", ""},
{Gnss::NormalEquationInfo::ESTIMATE_GRAVITY,                                "estimateGravity",                              "1", ""},
{Gnss::NormalEquationInfo::ESTIMATE_AMBIGUITIES,                            "estimateAmbiguities",                          "1", ""},
{Gnss::NormalEquationInfo::ESTIMATE_IONOSPHERE_STEC,                        "estimateIonosphereSTEC",                       "1", ""},
{Gnss::NormalEquationInfo::CONSTRAINT_IONOSPHERE_STEC,                      "constraintIonosphereSTEC",                     "1", ""},
{Gnss::NormalEquationInfo::ESTIMATE_IONOSPHERE_VTEC,                        "estimateIonosphereVTec",                       "1", ""},
{Gnss::NormalEquationInfo::ESTIMATE_IONOSPHERE_TECMAP,                      "estimateIonosphereTecMap",                     "1", ""},
{Gnss::NormalEquationInfo::CONSTRAINT_OTHER,                                "constraintOther",                              "1", ""}};

/***********************************************/

void Gnss::init(const std::vector<ReceiverPtr> &receiver, const std::vector<TransmitterPtr> &transmitter,
                const std::vector<ParametrizationPtr> &gnssParametrization,
                GnssParametrizationIonospherePtr gnssIonosphere, GnssParametrizationEarthRotationPtr gnssEarthRotation,
                GnssParametrizationGravityFieldPtr gnssGravityField, GnssParametrizationAmbiguitiesPtr gnssAmbiguities)
{
  try
  {
    this->receiver        = receiver;
    this->transmitter     = transmitter;
    this->parametrization = gnssParametrization;
    this->ionosphere      = gnssIonosphere;
    this->earthRotation   = gnssEarthRotation;
    this->gravityField    = gnssGravityField;
    this->ambiguities     = gnssAmbiguities;

    if(ionosphere)    this->parametrization.push_back(ionosphere);
    if(earthRotation) this->parametrization.push_back(earthRotation);
    if(gravityField)  this->parametrization.push_back(gravityField);
    if(ambiguities)   this->parametrization.push_back(ambiguities);

    for(UInt idRecv=0; idRecv<receiver.size(); idRecv++)
      receiver.at(idRecv)->registerGnss(this, idRecv);

    for(UInt idTrans=0; idTrans<transmitter.size(); idTrans++)
      transmitter.at(idTrans)->registerGnss(this, idTrans);

    for(UInt id=0; id<parametrization.size(); id++)
      parametrization.at(id)->_gnss = this;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::initInterval(AnalysisType analysisType, const std::vector<Time> &times, const Time &timeMargin, Parallel::CommunicatorPtr comm)
{
  try
  {
    this->times = times;

    for(const auto &para : parametrization) para->initIntervalEarly      (analysisType, times, timeMargin, comm);
    for(const auto &para : parametrization) para->initIntervalTransmitter(analysisType, times, timeMargin, comm);
    for(const auto &para : parametrization) para->initIntervalReceiver   (analysisType, times, timeMargin, comm);

    Bool disabled = FALSE;
    do
    {
      disabled = FALSE;
      for(const auto &para : parametrization)
        disabled = para->initIntervalDisabling(analysisType, times, timeMargin, comm) || disabled;
    }
    while(disabled);

    // distribute process id of station
    // --------------------------------
    Vector stationProcess(receiver.size());
    for(UInt idRecv=0; idRecv<receiver.size(); idRecv++)
      if(receiver.at(idRecv)->useable())
        stationProcess(idRecv) = Parallel::myRank(comm);
    Parallel::reduceSum(stationProcess, 0, comm);
    Parallel::broadCast(stationProcess, 0, comm);

    // collect observation types
    // -------------------------
    typesRecvTrans.clear();
    typesRecvTrans.resize(receiver.size(), std::vector<std::vector<GnssType>>(transmitter.size()));
    for(UInt idRecv=0; idRecv<receiver.size(); idRecv++)
    {
      if(receiver.at(idRecv)->useable())
        for(UInt idTrans=0; idTrans<transmitter.size(); idTrans++)
          typesRecvTrans.at(idRecv).at(idTrans) = receiver.at(idRecv)->observationsTypes(idTrans);
      Parallel::broadCast(typesRecvTrans.at(idRecv), static_cast<UInt>(stationProcess(idRecv)), comm);
    }

    // init other parametrizations
    // ---------------------------
    for(const auto &para : parametrization)
      para->initIntervalLate(analysisType, times, timeMargin, comm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::initParameter(NormalEquationInfo &normalEquationInfo)
{
  try
  {
    normalEquationInfo.initNewParameterNames();
    for(const auto &para : parametrization)
      para->initParameter(normalEquationInfo);
    normalEquationInfo.calculateIndex();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector Gnss::aprioriParameter(const NormalEquationInfo &normalEquationInfo) const
{
  try
  {
    Vector x0(normalEquationInfo.parameterCount());
    for(const auto &para : parametrization)
      para->aprioriParameter(normalEquationInfo, x0);
    Parallel::reduceSum(x0, 0, normalEquationInfo.comm);
    return x0;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool Gnss::isDesignMatrix(const NormalEquationInfo &normalEquationInfo, UInt idRecv, UInt idTrans, UInt idEpoch) const
{
  try
  {
    if(!normalEquationInfo.estimateReceiver.at(idRecv))
      return FALSE;

    // observation available?
    if(!receiver.at(idRecv)->observation(idTrans, idEpoch) || !transmitter.at(idTrans)->useable(idEpoch))
      return FALSE;

    if(receiver.at(idRecv)->isDesignMatrixReceiver(normalEquationInfo, idTrans, idEpoch))
      return TRUE;
    if(transmitter.at(idTrans)->isDesignMatrixTransmitter(normalEquationInfo, idRecv, idEpoch))
      return TRUE;
    for(const auto &para : parametrization)
      if(para->isDesignMatrix(normalEquationInfo, idRecv, idTrans, idEpoch))
        return TRUE;

    return FALSE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool Gnss::basicObservationEquations(const Gnss::NormalEquationInfo &normalEquationInfo, UInt idRecv, UInt idTrans, UInt idEpoch, Gnss::ObservationEquation &eqn) const
{
  try
  {
    if(!isDesignMatrix(normalEquationInfo, idRecv, idTrans, idEpoch))
      return FALSE;

    std::vector<GnssType> type;
    if(!receiver.at(idRecv)->observation(idTrans, idEpoch)->observationList(normalEquationInfo.analysisType, type))
      return FALSE;

    eqn.compute(*receiver.at(idRecv)->observation(idTrans, idEpoch), *receiver.at(idRecv), *transmitter.at(idTrans), ionosphere, idEpoch, TRUE, type);
    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::designMatrix(const NormalEquationInfo &normalEquationInfo, const ObservationEquation &eqn, DesignMatrix &A) const
{
  try
  {
    if(eqn.l.rows()==0)
      return;

    eqn.receiver->designMatrixReceiver(normalEquationInfo, eqn, A);
    eqn.transmitter->designMatrixTransmitter(normalEquationInfo, eqn, A);
    for(const auto &para : parametrization)
      para->designMatrix(normalEquationInfo, eqn, A);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::observationEquationEpoch(const NormalEquationInfo &normalEquationInfo, UInt idEpoch, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const
{
  try
  {
    for(const auto &para : parametrization)
      para->observationEquationEpoch(normalEquationInfo, idEpoch, normals, n, lPl, obsCount);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::observationEquation(const NormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const
{
  try
  {
    for(const auto &para : parametrization)
      para->observationEquation(normalEquationInfo, normals, n, lPl, obsCount);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double Gnss::ambiguityResolve(const NormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount, Bool dryRun)
{
  try
  {
    return ambiguities->ambiguityResolve(normalEquationInfo, normals, n, lPl, obsCount, dryRun);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double Gnss::updateParameter(const NormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz, Bool printStatistics)
{
  try
  {
    Double maxChange = 0;
    for(auto &para : parametrization)
      maxChange = std::max(maxChange, para->updateParameter(normalEquationInfo, x, Wz, printStatistics));
    return maxChange;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::updateCovariance(const NormalEquationInfo &normalEquationInfo, const MatrixDistributed &covariance)
{
  try
  {
    for(auto &para : parametrization)
      para->updateCovariance(normalEquationInfo, covariance);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::writeResults(const NormalEquationInfo &normalEquationInfo, const std::string &suffix)
{
  try
  {
    for(const auto &para : parametrization)
      para->writeResults(normalEquationInfo, suffix);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

std::vector<GnssType> Gnss::types(const GnssType mask) const
{
  try
  {
    std::vector<GnssType> types;
    for(UInt idRecv=0; idRecv<receiver.size(); idRecv++)
      for(UInt idTrans=0; idTrans<transmitter.size(); idTrans++)
        for(GnssType type : typesRecvTrans.at(idRecv).at(idTrans))
          if(GnssType::index(types, type) == NULLINDEX)
             types.push_back(type & mask);
    std::sort(types.begin(), types.end());
    return types;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::vector<GnssType> Gnss::replaceCompositeSignals(const std::vector<GnssType> &typesIn)
{
  try
  {
    std::vector<GnssType> types;

    for(GnssType t : typesIn)
      if((t == GnssType::PHASE) || (t == GnssType::RANGE)) // only phase and code signals are transmitted (what about doppler?)
      {
        GnssType prn = t & GnssType::PRN;
        if(t == GnssType::PHASE)
          types.push_back( t & ~GnssType::ATTRIBUTE );
        else if(t == GnssType::C2DG) {types.push_back(GnssType::C1CG + prn); types.push_back(GnssType::C1WG + prn); types.push_back(GnssType::C2WG + prn);}
        else if(t == GnssType::C1XG) {types.push_back(GnssType::C1SG + prn); types.push_back(GnssType::C1LG + prn);}
        else if(t == GnssType::C2XG) {types.push_back(GnssType::C2SG + prn); types.push_back(GnssType::C2LG + prn);}
        else if(t == GnssType::C5XG) {types.push_back(GnssType::C5IG + prn); types.push_back(GnssType::C5QG + prn);}
        else if(t == GnssType::C4XR) {types.push_back(GnssType::C4AR + prn); types.push_back(GnssType::C4BR + prn);}
        else if(t == GnssType::C6XR) {types.push_back(GnssType::C6AR + prn); types.push_back(GnssType::C6BR + prn);}
        else if(t == GnssType::C3XR) {types.push_back(GnssType::C3IR + prn); types.push_back(GnssType::C3QR + prn);}
  //       else if(t == GnssType::C1XE) {types.push_back(GnssType::C1BE + prn); types.push_back(GnssType::C1CE + prn);}
  //       else if(t == GnssType::C1ZE) {types.push_back(GnssType::C1AE + prn); types.push_back(GnssType::C1BE + prn); types.push_back(GnssType::C1CE + prn);}
  //       else if(t == GnssType::C5XE) {types.push_back(GnssType::C5IE + prn); types.push_back(GnssType::C5QE + prn);}
  //       else if(t == GnssType::C7XE) {types.push_back(GnssType::C7IE + prn); types.push_back(GnssType::C7QE + prn);}
  //       else if(t == GnssType::C8XE) {types.push_back(GnssType::C8IE + prn); types.push_back(GnssType::C8QE + prn);}
  //       else if(t == GnssType::C6XE) {types.push_back(GnssType::C6BE + prn); types.push_back(GnssType::C6CE + prn);}
  //       else if(t == GnssType::C6ZE) {types.push_back(GnssType::C6AE + prn); types.push_back(GnssType::C6BE + prn); types.push_back(GnssType::C6CE + prn);}
        else if(t == GnssType::C2XC) {types.push_back(GnssType::C2IC + prn); types.push_back(GnssType::C2QC + prn);}
        else if(t == GnssType::C1XC) {types.push_back(GnssType::C1DC + prn); types.push_back(GnssType::C1PC + prn);}
        else if(t == GnssType::C1ZC) {types.push_back(GnssType::C1SC + prn); types.push_back(GnssType::C1LC + prn);}
        else if(t == GnssType::C5XC) {types.push_back(GnssType::C5DC + prn); types.push_back(GnssType::C5PC + prn);}
        else if(t == GnssType::C7XC) {types.push_back(GnssType::C7IC + prn); types.push_back(GnssType::C7QC + prn);}
        else if(t == GnssType::C7ZC) {types.push_back(GnssType::C7DC + prn); types.push_back(GnssType::C7PC + prn);}
        else if(t == GnssType::C8XC) {types.push_back(GnssType::C8DC + prn); types.push_back(GnssType::C8PC + prn);}
        else if(t == GnssType::C6XC) {types.push_back(GnssType::C6IC + prn); types.push_back(GnssType::C6QC + prn);}
        else if(t == GnssType::C6ZC) {types.push_back(GnssType::C6DC + prn); types.push_back(GnssType::C6PC + prn);}
        // unknown attributes
        else if(t == GnssType::C2UG) {types.push_back(GnssType::C2SG + prn); types.push_back(GnssType::C2LG + prn);}
        else if(t == GnssType::C5UG) {types.push_back(GnssType::C5IG + prn); types.push_back(GnssType::C5QG + prn);}
        else if(t == GnssType::C1UE) {types.push_back(GnssType::C1XE + prn);}
        else if(t == GnssType::C5UE) {types.push_back(GnssType::C5XE + prn);}
        else if(t == GnssType::C7UE) {types.push_back(GnssType::C7XE + prn);}
        else if(t == GnssType::C8UE) {types.push_back(GnssType::C8XE + prn);}
        else if(t == GnssType::C6UE) {types.push_back(GnssType::C6XE + prn);}
  // ATTENTION: remove ? -> X replacement in defaultSignalComposition() if these composite types are reactivated
  //       else if(t == GnssType::C1UE) {types.push_back(GnssType::C1BE + prn); types.push_back(GnssType::C1CE + prn);}
  //       else if(t == GnssType::C5UE) {types.push_back(GnssType::C5IE + prn); types.push_back(GnssType::C5QE + prn);}
  //       else if(t == GnssType::C7UE) {types.push_back(GnssType::C7IE + prn); types.push_back(GnssType::C7QE + prn);}
  //       else if(t == GnssType::C8UE) {types.push_back(GnssType::C8IE + prn); types.push_back(GnssType::C8QE + prn);}
  //       else if(t == GnssType::C6UE) {types.push_back(GnssType::C6BE + prn); types.push_back(GnssType::C6CE + prn);}
        else
          types.push_back(t);
      }

    std::sort(types.begin(), types.end());
    types.erase(std::unique(types.begin(), types.end()), types.end()); // remove duplicates
    return types;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::defaultSignalComposition(const std::vector<GnssType> &types, std::vector<GnssType> &typesTrans, Matrix &A)
{
  try
  {
    // composed type = factor1 * type1 + factor2 * type2
    static const std::vector<std::tuple<GnssType, GnssType, GnssType, Double, Double>> composites =
      {{GnssType::C1XG,  GnssType::C1SG, GnssType::C1LG, 0.5, 0.5},
       {GnssType::C2XG,  GnssType::C2SG, GnssType::C2LG, 0.5, 0.5},
       {GnssType::C5XG,  GnssType::C5IG, GnssType::C5QG, 0.5, 0.5},
       {GnssType::C4XR,  GnssType::C4AR, GnssType::C4BR, 0.5, 0.5},
       {GnssType::C6XR,  GnssType::C6AR, GnssType::C6BR, 0.5, 0.5},
       {GnssType::C3XR,  GnssType::C3IR, GnssType::C3QR, 0.5, 0.5},
       {GnssType::C1XE,  GnssType::C1BE, GnssType::C1CE, 0.5, 0.5},
       {GnssType::C5XE,  GnssType::C5IE, GnssType::C5QE, 0.5, 0.5},
       {GnssType::C7XE,  GnssType::C7IE, GnssType::C7QE, 0.5, 0.5},
       {GnssType::C8XE,  GnssType::C8IE, GnssType::C8QE, 0.5, 0.5},
       {GnssType::C6XE,  GnssType::C6BE, GnssType::C6CE, 0.5, 0.5},

       {GnssType::C2XC,  GnssType::C2IC, GnssType::C2QC, 0.5, 0.5},
       {GnssType::C1XC,  GnssType::C1DC, GnssType::C1PC, 0.5, 0.5},
       {GnssType::C1ZC,  GnssType::C1SC, GnssType::C1LC, 0.5, 0.5},
       {GnssType::C5XC,  GnssType::C5DC, GnssType::C5PC, 0.5, 0.5},
       {GnssType::C7XC,  GnssType::C7IC, GnssType::C7QC, 0.5, 0.5},
       {GnssType::C7ZC,  GnssType::C7DC, GnssType::C7PC, 0.5, 0.5},
       {GnssType::C8XC,  GnssType::C8DC, GnssType::C8PC, 0.5, 0.5},
       {GnssType::C6XC,  GnssType::C6IC, GnssType::C6QC, 0.5, 0.5},
       {GnssType::C6ZC,  GnssType::C6DC, GnssType::C6PC, 0.5, 0.5},

       // unknown attributes
       {GnssType::C2UG,  GnssType::C2SG, GnssType::C2LG, 0.5, 0.5},
       {GnssType::C5UG,  GnssType::C5IG, GnssType::C5QG, 0.5, 0.5},
       {GnssType::C1UE,  GnssType::C1BE, GnssType::C1CE, 0.5, 0.5},
       {GnssType::C5UE,  GnssType::C5IE, GnssType::C5QE, 0.5, 0.5},
       {GnssType::C7UE,  GnssType::C7IE, GnssType::C7QE, 0.5, 0.5},
       {GnssType::C8UE,  GnssType::C8IE, GnssType::C8QE, 0.5, 0.5},
       {GnssType::C6UE,  GnssType::C6BE, GnssType::C6CE, 0.5, 0.5}};

    typesTrans = replaceCompositeSignals(types);

    A = Matrix(types.size(), typesTrans.size());
    for(UInt idType=0; idType<types.size(); idType++)
      if((types.at(idType) == GnssType::PHASE) || (types.at(idType) == GnssType::RANGE)) // only phase and code signals are transmitted (what about doppler?)
      {
        GnssType type = types.at(idType);

// ATTENTION: remove this if composite signals for Galileo are reactivated in replaceCompositeSignals()
if(type == GnssType::GALILEO + GnssType::RANGE + GnssType::UNKNOWN_ATTRIBUTE)
  type = (type & ~GnssType::ATTRIBUTE) + GnssType::X;

        const UInt idx = GnssType::index(typesTrans, type);
        if(idx != NULLINDEX)
        {
          A(idType, idx) = 1.; // signal observed directly
          continue;
        }

        if(type == GnssType::C2DG)
        {
          A(idType, GnssType::index(typesTrans, GnssType::C1CG)) = +1.;
          A(idType, GnssType::index(typesTrans, GnssType::C1WG)) = -1.;
          A(idType, GnssType::index(typesTrans, GnssType::C2WG)) = +1.;
          continue;
        }

        const auto composite = std::find_if(composites.begin(), composites.end(), [&](const auto &composite) {return (std::get<0>(composite) == type);});
        if(composite != composites.end())
        {
          A(idType, GnssType::index(typesTrans, std::get<1>(*composite))) = std::get<3>(*composite);
          A(idType, GnssType::index(typesTrans, std::get<2>(*composite))) = std::get<4>(*composite);
          continue;
        }

        throw(Exception("composite signal not implemented: "+type.str()));
      } // for(idType)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::vector<std::vector<GnssType>> Gnss::groupTypes(const std::vector<std::vector<std::vector<GnssType>>> &typesRecvTrans, UInt idDevice, Bool deviceIsReceiver) const
{
  try
  {
    const UInt transmitterSize = std::max_element(typesRecvTrans.begin(), typesRecvTrans.end(), [](const auto &a, const auto &b){return a.size() < b.size();})->size();
    std::vector<std::vector<GnssType>> groupTypes(deviceIsReceiver ? transmitterSize : typesRecvTrans.size());

    const UInt idRecvStart  = (deviceIsReceiver  && (idDevice != NULLINDEX)) ? idDevice : 0;
    const UInt idRecvEnd    = (deviceIsReceiver  && (idDevice != NULLINDEX)) ? idDevice : (typesRecvTrans.size()-1);
    const UInt idTransStart = (!deviceIsReceiver && (idDevice != NULLINDEX)) ? idDevice : 0;
    const UInt idTransEnd   = (!deviceIsReceiver && (idDevice != NULLINDEX)) ? idDevice : (transmitterSize-1);
    for(UInt idRecv=idRecvStart; idRecv<=idRecvEnd; idRecv++)
      for(UInt idTrans=idTransStart; idTrans<=idTransEnd; idTrans++)
        for(GnssType type : typesRecvTrans.at(idRecv).at(idTrans))
          if(GnssType::index(groupTypes.at(deviceIsReceiver ? idTrans : idRecv), type) == NULLINDEX)
            groupTypes.at(deviceIsReceiver ? idTrans : idRecv).push_back(type);

    for(UInt i=0; i<groupTypes.size(); i++)
      if(deviceIsReceiver)
        std::sort(groupTypes.at(i).begin(), groupTypes.at(i).end());
      else
        groupTypes.at(i) = replaceCompositeSignals(groupTypes.at(i));

    // delete empty groups
    groupTypes.erase(std::remove_if(groupTypes.begin(), groupTypes.end(), [](auto &t){return !t.size();}), groupTypes.end());

    // merge groups with identical range observations
    Bool merged;
    do
    {
      merged = FALSE;
      for(UInt i=0; i<groupTypes.size(); i++)
        for(UInt k=groupTypes.size()-1; k>i; k--)
        {
          Bool identical = FALSE;
          for(GnssType type : groupTypes.at(k))
            if((type == GnssType::RANGE) && (GnssType::index(groupTypes.at(i), type) != NULLINDEX))
            {
              identical = TRUE;
              break;
            }

          if(identical)
          {
            for(GnssType type : groupTypes.at(k))
              if(GnssType::index(groupTypes.at(i), type) == NULLINDEX)
                groupTypes.at(i).push_back(type);
            std::sort(groupTypes.at(i).begin(), groupTypes.at(i).end());
            groupTypes.erase(groupTypes.begin()+k);
            merged = TRUE;
          }
        } // for(group i,k)
    }
    while(merged);

    return groupTypes;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::signalBiasParameter(const std::vector<GnssType> &obsTypes, const std::vector<std::vector<GnssType>> &groupTypes,
                               Bool eliminateClock, Bool estimateCodeBias, Bool estimatePhaseBias, Bool estimateTecBias,
                               Matrix &A, std::vector<ParameterName> &parameterName) const
{
  try
  {
    // --- lambda ---------------------
    auto extendColumns = [](Matrix &A, UInt columns)
    {
      Matrix tmp = A;
      A = Matrix(tmp.rows(), tmp.columns()+columns);
      copy(tmp, A.column(0, tmp.columns()));
      return tmp.columns();
    };

    A = Matrix(obsTypes.size(), 0);

    // Code biases
    // -----------
    if(estimateCodeBias)
    {
      std::vector<GnssType> biasTypes;
      for(auto &group : groupTypes)
        for(auto type : group)
          if((type == GnssType::RANGE) && (GnssType::index(obsTypes, type) != NULLINDEX) && (GnssType::index(biasTypes, type) == NULLINDEX))
            biasTypes.push_back(type);

      // transformation matrix
      A = Matrix(obsTypes.size(), biasTypes.size());
      for(UInt idType=0; idType<biasTypes.size(); idType++)
        A(GnssType::index(obsTypes, biasTypes.at(idType)), idType) = 1.;

      // determine null space of clock per system and TEC
      Matrix B(biasTypes.size(), 1+eliminateClock*groupTypes.size());
      const UInt idxC1WG = GnssType::index(biasTypes, GnssType::C1WG);
      const UInt idxC2WG = GnssType::index(biasTypes, GnssType::C2WG);
      if((groupTypes.size() == 1) && (idxC1WG != NULLINDEX) && (idxC2WG != NULLINDEX))  // special case: GPS clock/TEC is defined via C1WG/C2WG only
      {
        B(idxC1WG, 0) = GnssType::C1WG.ionosphericFactor();
        B(idxC2WG, 0) = GnssType::C2WG.ionosphericFactor();
        if(eliminateClock)
          B(idxC1WG, 1) = B(idxC2WG, 1) = 1.;
      }
      else
        for(UInt idGroup=0; idGroup<groupTypes.size(); idGroup++)
          for(auto type : groupTypes.at(idGroup))
          {
            const UInt idx = GnssType::index(biasTypes, type);
            if(idx == NULLINDEX)
              continue;
            B(idx, 0) = type.ionosphericFactor();
            if(eliminateClock)
              B(idx, 1+idGroup) = 1;
          }

      // eliminate nullspace (clock and TECs)
      const Vector tau = QR_decomposition(B);
      Matrix I = identityMatrix(biasTypes.size());
      QMult(B, tau, I);
      A = A * I.column(B.columns(), I.columns()-B.columns());
    } // if(estimateCodeBias)

    // Phase biases
    // ------------
    if(estimatePhaseBias)
    {
      std::vector<GnssType> biasTypes;
      for(auto &group : groupTypes)
        for(auto type : group)
          if((type == GnssType::PHASE) && (GnssType::index(obsTypes, type) != NULLINDEX) && (GnssType::index(biasTypes, type) == NULLINDEX))
            biasTypes.push_back(type);
      std::sort(biasTypes.begin(), biasTypes.end());
      // transformation matrix
      const UInt column = extendColumns(A, biasTypes.size());
      for(UInt idType=0; idType<biasTypes.size(); idType++)
        A(GnssType::index(obsTypes, biasTypes.at(idType)), column + idType) = 1.;
    }

    // determine parameter names
    // -------------------------
    for(UInt i=0; i<A.columns(); i++)
    {
      std::string typeStr;
      for(UInt idType=0; idType<A.rows(); idType++)
        if(std::fabs(A(idType, i)) > 1e-4)
          typeStr += ((A(idType, i) > 0) ? "+" : "") + A(idType, i)%"%.2f"s + obsTypes.at(idType).str();
      parameterName.push_back(ParameterName("", "signalBias"+(i+1)%"%02i("s+typeStr+")"));
    }

    // TEC bias parameter
    // ------------------
    if(estimateTecBias)
    {
      // test for common types
      Bool common = FALSE;
      for(GnssType type : obsTypes)
        common = common || (std::count_if(groupTypes.begin(), groupTypes.end(), [&](const auto &types){return (GnssType::index(types, type) != NULLINDEX);}) > 1);

      UInt countTec = 1;
      if(!common)
        countTec = groupTypes.size();

      // transformation matrix
      const UInt column = extendColumns(A, countTec);
      for(UInt idType=0; idType<obsTypes.size(); idType++)
        for(UInt idGroup=0; idGroup<groupTypes.size(); idGroup++)
          if(GnssType::index(groupTypes.at(idGroup), obsTypes.at(idType)) != NULLINDEX)
            A(idType, column+std::min(idGroup, countTec-1)) = obsTypes.at(idType).ionosphericFactor();

      for(UInt i=0; i<countTec; i++)
        parameterName.push_back(ParameterName("", "signalBias"+(i+1)%"%02iTEC"s));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::checkMaxChange(Double &maxChange, std::string &info, Bool printStatistics, Parallel::CommunicatorPtr comm)
{
  try
  {
    Vector change(Parallel::size(comm));
    change(Parallel::myRank(comm)) = maxChange;
    Parallel::reduceSum(change, 0, comm);
    Parallel::broadCast(change, 0, comm);
    UInt idProcess = 0;
    maxChange = 0;
    for(UInt i=0; i<change.rows(); i++)
      if(change(i) > maxChange)
      {
        maxChange = change(i);
        idProcess = i;
      }

    if(printStatistics && (idProcess != 0) && (idProcess == Parallel::myRank(comm)))
      Parallel::send(info, 0, comm);
    else if(printStatistics && Parallel::isMaster(comm))
    {
      if(idProcess != 0)
        Parallel::receive(info, idProcess, comm);
      if(!info.empty())
        logInfo<<info<<Log::endl;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

Gnss::NormalEquationInfo::NormalEquationInfo(UInt countEpoch, UInt countReceiver, UInt countTransmitter, Parallel::CommunicatorPtr comm_) :
    comm(comm_),
    analysisType(ANALYSIS_RAW),
    estimationType(MASK_ALL),
    estimateReceiver(countReceiver, TRUE),
    idEpochs(countEpoch),
    defaultBlockSizeEpoch(0),
    defaultBlockSizeInterval(64),
    defaultBlockSizeAmbiguity(64),
    defaultBlockReceiverCount(0),
    defaultBlockCountReduction(32),
    keepEpochNormalsInMemory(FALSE),
    accumulateEpochObservations(FALSE),
    blockCountEpoch_(countEpoch, 0),
    countTransmitter_(countTransmitter)
{
  std::iota(idEpochs.begin(), idEpochs.end(), 0);
}

/***********************************************/

void Gnss::NormalEquationInfo::initNewParameterNames()
{
  try
  {
    parameters.clear();
    block_.clear();
    index_.clear();
    count_.clear();
    parameterNames_.clear();
    blockIndices_.clear();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Gnss::ParameterIndex Gnss::NormalEquationInfo::addParameters(UInt idEpoch, UInt idRecv, UInt idTrans, const std::vector<ParameterName> &parameterNames)
{
  if(!parameterNames.size())
    return ParameterIndex(NULLINDEX);
  parameters.push_back(std::make_tuple(idEpoch, idRecv, idTrans, parameters.size(), parameterNames));
  return ParameterIndex(parameters.size()-1);
}

/***********************************************/

void Gnss::NormalEquationInfo::calculateIndex()
{
  try
  {
    auto newBlock = [&]()
    {
      const UInt parameterCount = this->parameterCount();
      if(blockSize(blockCount()-1))
        blockIndices_.push_back(parameterCount);
    };

    auto insert = [&](auto iter, UInt defaultBlockSize)
    {
      if((defaultBlockSize > 0) && (blockSize(blockCount()-1) >= defaultBlockSize))
        newBlock();
      const UInt idx   = std::get<3>(*iter);
      index_.at(idx)   = blockIndices_.back();
      block_.at(idx)   = blockCount()-1;
      count_.at(idx)   = std::get<4>(*iter).size();
      parameterNames_.insert(parameterNames_.end(), std::get<4>(*iter).begin(), std::get<4>(*iter).end());
      blockIndices_.back() += std::get<4>(*iter).size();
    };

    block_.resize(parameters.size(), NULLINDEX);
    index_.resize(parameters.size(), NULLINDEX);
    count_.resize(parameters.size(), 0);
    blockIndices_ = {0, 0};
    parameterNames_.reserve(std::accumulate(parameters.begin(), parameters.end(), UInt(0), [](UInt count, const auto &p){return count+std::get<4>(p).size();}));

    std::stable_sort(parameters.begin(), parameters.end(), [](auto &p1, auto &p2)
                    {
                      if((std::get<1>(p1) != NULLINDEX) && (std::get<2>(p1) != NULLINDEX)) return FALSE; // ambiguities always at end
                      if(std::get<0>(p1) != std::get<0>(p2)) return (std::get<0>(p1) < std::get<0>(p2)); // epoch
                      if(std::get<1>(p1) != std::get<1>(p2)) return (std::get<1>(p1) < std::get<1>(p2)); // idRecv
                      return (std::get<2>(p1) < std::get<2>(p2));                                        // idTrans
                    });
    auto iter = parameters.begin();

    // epoch parameters
    std::fill(blockCountEpoch_.begin(), blockCountEpoch_.end(), 0);
    for(UInt idEpoch : idEpochs)
    {
      newBlock();
      UInt blockEpochStart = blockIndices_.size();
      UInt idRecv = NULLINDEX;
      UInt countStation = 0;
      while((iter != parameters.end()) && (std::get<0>(*iter) == idEpoch) && ((std::get<1>(*iter) == NULLINDEX) || (std::get<2>(*iter) == NULLINDEX)))
      {
        if(std::get<1>(*iter) != idRecv) // next receiver?
          if(defaultBlockReceiverCount && ((std::get<1>(*iter) == NULLINDEX) || ((countStation++ % defaultBlockReceiverCount) == 0)))
            newBlock();
        idRecv = std::get<1>(*iter);

        insert(iter, defaultBlockSizeEpoch);
        iter++;
      }
      blockCountEpoch_.at(idEpoch) = blockIndices_.size() - blockEpochStart + (blockSize(blockCount()-1) ? 1 : 0);
    }

    // receiver interval parameters
    newBlock();
    blockInterval_ = blockCount()-1;
    UInt countStation = 0;
    for(UInt idRecv=0; idRecv<estimateReceiver.size(); idRecv++)
    {
      if(defaultBlockReceiverCount && ((countStation++ % defaultBlockReceiverCount) == 0))
        newBlock();
      Bool firstBlock = TRUE;
      while((iter != parameters.end()) && (std::get<0>(*iter) == NULLINDEX) && (std::get<1>(*iter) == idRecv) && (std::get<2>(*iter) == NULLINDEX))
      {
        insert(iter, firstBlock ? defaultBlockSizeInterval : 0); // do not split parameters of a receiver
        firstBlock = FALSE;
        iter++;
      }
    }

    // transmitter interval parameters
    newBlock();
    for(UInt idTrans=0; idTrans<countTransmitter_; idTrans++)
    {
      Bool firstBlock = TRUE;
      while((iter != parameters.end()) && (std::get<0>(*iter) == NULLINDEX) && (std::get<1>(*iter) == NULLINDEX) && (std::get<2>(*iter) == idTrans))
      {
        insert(iter, firstBlock ? defaultBlockSizeInterval : 0); // do not split parameters of a transmitter
        firstBlock = FALSE;
        iter++;
      }
    }

    // other interval parameters
    newBlock();
    while((iter != parameters.end()) && (std::get<0>(*iter) == NULLINDEX) && (std::get<1>(*iter) == NULLINDEX) && (std::get<2>(*iter) == NULLINDEX))
    {
      insert(iter, defaultBlockSizeInterval);
      iter++;
    }

    // ambiguity parameters
    newBlock();
    blockAmbiguity_ = blockCount()-1;
    while((iter != parameters.end()) && (std::get<1>(*iter) != NULLINDEX) && (std::get<2>(*iter) != NULLINDEX))
    {
      insert(iter, defaultBlockSizeAmbiguity);
      iter++;
    }

    // remove possible last empty block
    if(blockCount() && !blockSize(blockCount()-1))
      blockIndices_.pop_back();

    parameters.clear();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::NormalEquationInfo::removeEpochParameters()
{
  try
  {
    parameterNames_.erase(parameterNames_.begin(), parameterNames_.begin()+blockInterval());
    blockIndices_.erase(blockIndices_.begin(), blockIndices_.begin()+blockInterval());
    for(UInt &blockIndex : blockIndices_)
      blockIndex -= blockInterval();
    blockAmbiguity_ -= blockInterval();
    blockInterval_   = 0;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

Gnss &Gnss::Parametrization::gnss() const
{
  if(_gnss == nullptr)
    throw(Exception("Parametrization is not registered in Gnss class"));
  return *_gnss;
}

/***********************************************/
