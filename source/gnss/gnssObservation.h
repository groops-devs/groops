/***********************************************/
/**
* @file gnssObservation.h
*
* @brief Code & Phase observations.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2012-04-18
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSOBSERVATION__
#define __GROOPS_GNSSOBSERVATION__

#include "base/gnssType.h"

/** @addtogroup gnssGroup */
/// @{

/***** TYPES ***********************************/

class GnssTrack;
class GnssReceiver;
class GnssTransmitter;

/***** CLASS ***********************************/

/** @brief Single observation. */
class GnssSingleObservation
{
public:
  GnssType type;        ///< measurement types (phases, pseudo ranges, ...)
  Double   observation; ///< original observations
  Double   residuals;   ///< estimated postfit residuals
  Double   redundancy;  ///< partial redundancies of the least squares adjustment
  Double   sigma0;      ///< expected (apriori) accuracies
  Double   sigma;       ///< modified accuracies (downweighted outliers)

  GnssSingleObservation() {}
  GnssSingleObservation(GnssType _type, Double _observation=0., Double _residuals=0., Double _redundancy=0., Double _sigma0=0., Double _sigma=0.)
    : type(_type), observation(_observation), residuals(_residuals), redundancy(_redundancy), sigma0(_sigma0), sigma(_sigma) {}
};

/***** CLASS ***********************************/

/** @brief Observations.
* Between one receiver and one transmitter at one epoch. */
class GnssObservation
{
  std::vector<GnssSingleObservation> obs;

public:
  typedef UInt Group;
  constexpr static Group RANGE = 1<<0; // use code observations only
  constexpr static Group PHASE = 1<<1;
  constexpr static Group IONO  = 1<<2;
  constexpr static Group ALL   = ~0;   // all bits set: use all available observations as they are

  GnssTrack *track; ///< phase ambiguities
  Double     STEC, dSTEC;  ///< total ionosphere slant TEC along the path and estimated part of STEC

  GnssObservation() : track(nullptr), STEC(0.), dSTEC(0.) {}

  UInt size() const                                      {return obs.size();}
  GnssSingleObservation       &at(UInt idType)           {return obs.at(idType);}
  const GnssSingleObservation &at(UInt idType) const     {return obs.at(idType);}
  GnssSingleObservation       &at(GnssType type)         {return obs.at(index(type));}
  const GnssSingleObservation &at(GnssType type) const   {return obs.at(index(type));}
  UInt index(GnssType type) const;
  void push_back(const GnssSingleObservation &singleObs) {obs.push_back(singleObs);}
  void erase(UInt idType)                                {obs.erase(obs.begin()+idType); obs.shrink_to_fit();}
  void shrink_to_fit()                                   {obs.shrink_to_fit();}

  Bool init(const GnssReceiver &receiver, const GnssTransmitter &transmitter, const std::function<Rotary3d(const Time &time)> &rotationCrf2Trf,
            UInt idEpoch, Angle elevationCutOff, Double &phaseWindupOld);

  /** @brief Returns true if observation types required by @a analysisType are available, @a type contains list of these observations. */
  Bool observationList         (Group group, std::vector<GnssType> &types) const;
  void setDecorrelatedResiduals(const std::vector<GnssType> &types, const_MatrixSliceRef residuals, const_MatrixSliceRef redundancy);
  void updateParameter         (const_MatrixSliceRef x, const_MatrixSliceRef covariance=Matrix());
};

/***** CLASS ***********************************/

/** @brief Reduced observations (obs - computed) and design matrix.
* Between one receiver and one transmitter at one epoch. */
class GnssObservationEquation
{
public:
  enum {idxPosRecv    = 0, // x,y,z (CRF)
        idxClockRecv  = 3,
        idxPosTrans   = 4, // x,y,z (CRF)
        idxClockTrans = 7,
        idxRange      = 8,
        idxSTEC       = 9,
        idxUnit       = 10};

  UInt  idEpoch;
  const GnssTrack       *track; // phase ambiguities
  const GnssReceiver    *receiver;
  const GnssTransmitter *transmitter;

  // weighted observations (with 1/sigma)
  std::vector<GnssType> types;            ///< observed types (inclusive composed signals)
  std::vector<GnssType> typesTransmitted; ///< original transmitted signals (C2XG -> C2LG + C2SG), phases without attribute
  UInt   rankDeficit;  ///< from eliminated group parameters
  Vector l;            ///< weighted reduced observations
  Vector sigma;
  Vector sigma0;

  // design matrix
  Matrix A;      ///< columns: dl/dx, dl/dy, dl/dz, dl/dClock, unit matrix, transformation matrix (typesTransmitter->types)
  Matrix B;      ///< ionosphere, ...

  // approximate values (Taylor point)
  Time     timeRecv, timeTrans;
  Vector3d posRecv,  posTrans;
  Vector3d velocityRecv, velocityTrans;
  Angle    azimutRecvLocal, elevationRecvLocal;
  Angle    azimutRecvAnt,   elevationRecvAnt;
  Angle    azimutTrans,     elevationTrans;
  Double   STEC, dSTEC;

  GnssObservationEquation() : idEpoch(NULLINDEX), track(nullptr), receiver(nullptr), transmitter(nullptr), rankDeficit(0), STEC(0), dSTEC(0) {}

  GnssObservationEquation(const GnssObservation &observation, const GnssReceiver &receiver, const GnssTransmitter &transmitter,
                          const std::function<Rotary3d(const Time &time)> &rotationCrf2Trf, const std::function<void(GnssObservationEquation &eqn)> &reduceModels,
                          UInt idEpoch, Bool decorrelate, const std::vector<GnssType> &types)
    {compute(observation, receiver, transmitter, rotationCrf2Trf, reduceModels, idEpoch, decorrelate, types);}

  void compute(const GnssObservation &observation, const GnssReceiver &receiver, const GnssTransmitter &transmitter,
               const std::function<Rotary3d(const Time &time)> &rotationCrf2Trf, const std::function<void(GnssObservationEquation &eqn)> &reduceModels,
               UInt idEpoch, Bool decorrelate, const std::vector<GnssType> &types);

  void eliminateGroupParameters(Bool removeRows=TRUE);
};

/***********************************************/

/// @}

#endif /* __GROOPS___ */
