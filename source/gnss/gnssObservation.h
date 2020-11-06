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

#include "gnss/gnss.h"

/** @addtogroup gnssGroup */
/// @{

/***** CLASS ***********************************/

/** @brief Single observation. */
class Gnss::SingleObservation
{
public:
  GnssType type;        ///< measurement types (phases, pseudo ranges, ...)
  Double   observation; ///< original observations
  Double   residuals;   ///< estimated postfit residuals
  Double   redundancy;  ///< partial redundancies of the least squares adjustment
  Double   sigma0;      ///< expected (apriori) accuracies
  Double   sigma;       ///< modfied accuracies (downweighted outliers)

  SingleObservation() {}
  SingleObservation(GnssType _type, Double _observation=0., Double _residuals=0., Double _redundancy=0., Double _sigma0=0., Double _sigma=0.)
  : type(_type), observation(_observation), residuals(_residuals), redundancy(_redundancy), sigma0(_sigma0), sigma(_sigma) {}
};

/***** CLASS ***********************************/

/** @brief Observations.
* Between one receiver and one transmitter at one epoch. */
class Gnss::Observation
{
  std::vector<SingleObservation> obs;

public:
  Track  *track;     ///< phase ambiguities
  Double  dSTEC;     ///< estimated ionosphere slant TEC along the path (additional to model)

  Observation() : track(nullptr), dSTEC(0.) {}

  UInt size() const                                  {return obs.size();}
  SingleObservation       &at(UInt idType)           {return obs.at(idType);}
  const SingleObservation &at(UInt idType) const     {return obs.at(idType);}
  SingleObservation       &at(GnssType type)         {return obs.at(index(type));}
  const SingleObservation &at(GnssType type) const   {return obs.at(index(type));}
  void push_back(const SingleObservation &singleObs) {obs.push_back(singleObs);}
  void erase(UInt idType)                            {obs.erase(obs.begin()+idType); obs.shrink_to_fit();}
  void shrink_to_fit()                               {obs.shrink_to_fit();}

  Bool init(const Receiver &receiver, const Transmitter &transmitter, GnssParametrizationIonospherePtr ionosphere, UInt idEpoch, Double &phaseWindupOld);
  UInt index(GnssType type) const;

  /** @brief Returns true if observation types required by @a analysisType are available, @a type contains list of these observations. */
  Bool observationList         (AnalysisType analysisType, std::vector<GnssType> &types) const;
  void setDecorrelatedResiduals(const std::vector<GnssType> &types, const_MatrixSliceRef residuals, const_MatrixSliceRef redundancy);
  void updateParameter         (const_MatrixSliceRef x, const_MatrixSliceRef covariance=Matrix());
};

/***** CLASS ***********************************/

/** @brief Reduced observations (obs - computed) and design matrix.
* Between one receiver and one transmitter at one epoch. */
class Gnss::ObservationEquation
{
public:
  enum {idxPosRecv    = 0, // x,y,z (CRF)
        idxClockRecv  = 3,
        idxPosTrans   = 4, // x,y,z (CRF)
        idxClockTrans = 7,
        idxRange      = 8,
        idxUnit       = 9};

  UInt  idEpoch;
  const Track       *track; // phase ambiguities
  const Receiver    *receiver;
  const Transmitter *transmitter;

  // weighted observations (with 1/sigma)
  std::vector<GnssType> types;              // observed types (inclusive composed signals)
  std::vector<GnssType> typesTransmitted;   // orginal transmitted signals (C2XG -> C2LG + C2SG), phases without attribute
  Vector l;      ///< weighted reduced observations
  Vector sigma;
  Vector sigma0;

  // design matrix
  Matrix A;      ///< columns: dl/dx, dl/dy, dl/dz, dl/dClock, unit matrix, transformation matrix (typesTransmitter->types)
  Matrix B;      ///< ionosphere, ...

  // approximate values (Taylor point)
  Time     timeRecv, timeTrans;
  Vector3d posRecv, posTrans;
  Vector3d velocityRecv, velocityTrans;
  Double   r12; // corrected range: range - (clockError, troposphere, relativistic effects)
  Angle    azimutRecvLocal, elevationRecvLocal;
  Angle    azimutRecvAnt,   elevationRecvAnt;
  Angle    azimutTrans,     elevationTrans;
  Double   dSTEC;

  ObservationEquation() : track(nullptr), receiver(nullptr), transmitter(nullptr) {}

  ObservationEquation(const Observation &observation, const Receiver &receiver, const Transmitter &transmitter,
                      GnssParametrizationIonospherePtr ionosphere, UInt idEpoch, Bool decorrelate, const std::vector<GnssType> &types)
                      {compute(observation, receiver, transmitter, ionosphere, idEpoch, decorrelate, types);}

  void compute(const Observation &observation, const Receiver &receiver, const Transmitter &transmitter,
               GnssParametrizationIonospherePtr ionosphere, UInt idEpoch, Bool decorrelate, const std::vector<GnssType> &types);
};

/***********************************************/

/// @}

#endif /* __GROOPS___ */
