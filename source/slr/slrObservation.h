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

  void setHomogenizedResiduals(const_MatrixSliceRef residuals, const_MatrixSliceRef redundancy);
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

  enum Type {RANGE, DIRECTIONS};

  Type  type;
  const SlrStation   *station;
  const SlrSatellite *satellite;

  // weighted observations (with 1/sigma)
  Vector l;
  Matrix A; // design matrix
  Vector sigmas;
  Vector sigmas0;

  // approximate values (Taylor point)
  std::vector<UInt>     index, count; // observations per epoch
  std::vector<Vector3d> posStat,   posSat;
  std::vector<Time>     timesStat, timesSat;
  std::vector<Angle>    azimutStat, elevationStat;
  Vector                laserWavelength;

  SlrObservationEquation() : station(nullptr), satellite(nullptr) {}

  void compute(const SlrObservation &observation, const SlrStation &station, const SlrSatellite &satellite,
               const std::function<Rotary3d(const Time &time)> &rotationCrf2Trf,
               const std::function<void(SlrObservationEquation &eqn)> &reduceModels, Bool homogenize);
};

/***********************************************/

/// @}

#endif /* __GROOPS___ */
