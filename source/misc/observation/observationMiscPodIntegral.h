/***********************************************/
/**
* @file observationMiscPodIntegral.h
*
* @brief Precise Orbit data (Short Arc Integral).
*
* @author Torsten Mayer-Guerr
* @date 2015-06-02
*
*/
/***********************************************/

#ifndef __GROOPS_OBSERVATIONMSICPODINTEGRAL__
#define __GROOPS_OBSERVATIONMSICPODINTEGRAL__

/***********************************************/

#include "files/fileInstrument.h"
#include "files/fileSatelliteModel.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/ephemerides/ephemerides.h"
#include "classes/gravityfield/gravityfield.h"
#include "classes/parametrizationGravity/parametrizationGravity.h"
#include "classes/parametrizationAcceleration/parametrizationAcceleration.h"
#include "misc/observation/integralEquation.h"
#include "misc/observation/observationMisc.h"
#include "misc/observation/observationMiscPod.h"

/***** TYPES ***********************************/

class ObservationMiscPodIntegral;
typedef std::shared_ptr<ObservationMiscPodIntegral> ObservationMiscPodIntegralPtr;

/***** CLASS ***********************************/

/** @brief Precise Orbit data (Short Arc Integral).
* @ingroup miscGroup
* @see Observation */
class ObservationMiscPodIntegral : public ObservationMiscPod
{
  SatelliteModelPtr              satellite;
  std::vector<PodRightSidePtr>   rhs; // right hand sides
  IntegralEquation               integralEquation;
  Bool                           keepSatelliteStates;
  Bool                           accelerateComputation;
  InstrumentFile                 orbitFile;
  InstrumentFile                 starCameraFile;
  EarthRotationPtr               earthRotation;
  EphemeridesPtr                 ephemerides;
  GravityfieldPtr                gradientfield;
  ParametrizationGravityPtr      parameterGravity;
  ParametrizationAccelerationPtr parameterAcceleration;
  UInt                           integrationDegree;
  UInt                           interpolationDegree;
  UInt                           countArc;

  // Indicies for design matrix A
  UInt countAParameter;
  UInt idxGravity;
  UInt idxBound;
  UInt idxSat;

public:
  ObservationMiscPodIntegral(Config &config);
 ~ObservationMiscPodIntegral() {}

  Bool setInterval(const Time &timeStart, const Time &timeEnd) override;
  UInt parameterCount()          const override {return countAParameter;}
  UInt gravityParameterCount()   const override {return parameterGravity->parameterCount();}
  UInt rightSideCount()          const override {return rhs.size();}
  UInt arcCount()                const override {return countArc;}
  void parameterName(std::vector<ParameterName> &name) const override;

  Arc computeArc(UInt arcNo, CovariancePodPtr covPod);
};

/***********************************************/

#endif /* __GROOPS__ */
