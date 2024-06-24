/***********************************************/
/**
* @file slrObservation.h
*
* @brief Code & Phase observations.
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#ifndef __GROOPS_SLROBSERVATION__
#define __GROOPS_SLROBSERVATION__

/** @addtogroup slrGroup */
/// @{

/***** TYPES ***********************************/

class SlrStation;
class SlrSatellite;
class SlrObservation;
typedef std::shared_ptr<SlrObservation> SlrObservationPtr;

/***** CLASS ***********************************/

/** @brief Observations.
* Between one station and one satellite (one pass). */
class SlrObservation
{
public:
  std::vector<Time> timesTrans;      ///< transmitted time
  Vector            observations;    ///< original observations
  Vector            residuals;       ///< estimated postfit residuals
  Vector            redundancies;    ///< partial redundancies of the least squares adjustment
  Vector            sigmas0;         ///< expected (apriori) accuracies
  Vector            sigmas;          ///< modfied accuracies (downweighted outliers)
  Vector            laserWavelength; ///< laser wavelength

  SlrObservation() {}

  Bool init(const SlrStation &station, const SlrSatellite &satellite,
            const std::function<Rotary3d(const Time &time)> &rotationCrf2Trf, Angle elevationCutOff);

  void setDecorrelatedResiduals(const_MatrixSliceRef residuals, const_MatrixSliceRef redundancy);
};

/***** CLASS ***********************************/

/** @brief Reduced observations (obs - computed) and design matrix.
* Between one station and one satellite for one pass. */
class SlrObservationEquation
{
public:
  enum {idxPosStat = 0,  // x,y,z (CRF)
        idxPosSat  = 3,  // x,y,z (CRF)
        idxTime    = 6,  // transmit time offset
        idxRange   = 7}; // One way range

  const SlrObservation *observation;
  const SlrStation     *station;
  const SlrSatellite   *satellite;

  // weighted observations (with 1/sigma)
  Vector l;
  Matrix A; // design matrix
  Vector sigmas;
  Vector sigmas0;

  // approximate values (Taylor point)
  std::vector<Time>     timesTrans, timesBounce, timesRecv;
  std::vector<Vector3d> posTrans, posBounce, posRecv;
  std::vector<Vector3d> velTrans, velBounce;
  std::vector<Angle>    azimutStat, elevationStat;
  Vector                laserWavelength;

  SlrObservationEquation() : observation(nullptr), station(nullptr), satellite(nullptr) {}

  SlrObservationEquation(const SlrObservation &observation, const SlrStation &station, const SlrSatellite &satellite,
                         const std::function<Rotary3d(const Time &time)> &rotationCrf2Trf,
                         const std::function<void(SlrObservationEquation &eqn)> &reduceModels, Bool decorrelate)
    {compute(observation, station, satellite, rotationCrf2Trf, reduceModels, decorrelate);}

  void compute(const SlrObservation &observation, const SlrStation &station, const SlrSatellite &satellite,
               const std::function<Rotary3d(const Time &time)> &rotationCrf2Trf,
               const std::function<void(SlrObservationEquation &eqn)> &reduceModels, Bool decorrelate);
};

/***********************************************/

/// @}

#endif /* __GROOPS___ */
