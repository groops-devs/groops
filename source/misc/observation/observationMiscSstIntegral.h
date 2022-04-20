/***********************************************/
/**
* @file observationMiscSstIntegral.h
*
* @brief Satellite to satellite tracking (Short Arc Integral).
*
* @author Torsten Mayer-Guerr
* @date 2015-06-02
*
*/
/***********************************************/

#ifndef __GROOPS_OBSERVATIONMSICSSTINTEGRAL__
#define __GROOPS_OBSERVATIONMSICSSTINTEGRAL__

/***********************************************/

#include "files/fileInstrument.h"
#include "files/fileSatelliteModel.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/ephemerides/ephemerides.h"
#include "classes/gravityfield/gravityfield.h"
#include "classes/parametrizationGravity/parametrizationGravity.h"
#include "classes/parametrizationAcceleration/parametrizationAcceleration.h"
#include "classes/parametrizationSatelliteTracking/parametrizationSatelliteTracking.h"
#include "misc/observation/integralEquation.h"
#include "misc/observation/covarianceSst.h"
#include "misc/observation/covariancePod.h"
#include "misc/observation/observationMisc.h"
#include "misc/observation/observationMiscSst.h"

/***** TYPES ***********************************/

class ObservationMiscSstIntegral;
typedef std::shared_ptr<ObservationMiscSstIntegral> ObservationMiscSstIntegralPtr;

/***** CLASS ***********************************/

/** @brief Satellite to satellite tracking (Short Arc Integral).
* @ingroup miscGroup
* @see Observation */
class ObservationMiscSstIntegral : public ObservationMiscSst
{
  SatelliteModelPtr                   satellite1;
  SatelliteModelPtr                   satellite2;
  std::vector<SstRightSidePtr>        rhs; // right hand sides
  IntegralEquation                    integralEquation;
  Bool                                keepSatelliteStates;
  InstrumentFile                      orbit1File, orbit2File;
  InstrumentFile                      starCamera1File, starCamera2File;
  EarthRotationPtr                    earthRotation;
  EphemeridesPtr                      ephemerides;
  GravityfieldPtr                     gradientfield;
  ParametrizationGravityPtr           parameterGravity;
  ParametrizationAccelerationPtr      parameterAcceleration1;
  ParametrizationAccelerationPtr      parameterAcceleration2;
  ParametrizationSatelliteTrackingPtr parameterSst;
  UInt                                integrationDegree;
  UInt                                interpolationDegree;
  UInt                                countArc;
  UInt                                sstType; // 0: biased range, 1: range-rate, 2: range-acceleration
  Bool                                computeRange, computeVelocity, computeAcceleration;

  // Indicies for design matrix A
  UInt countAParameter;
  UInt idxGravity, gravityCount;
  UInt idxBound1,  idxBound2;
  UInt idxSat1, idxSat2;
  UInt idxSstPara;

public:
  ObservationMiscSstIntegral(Config &config);
 ~ObservationMiscSstIntegral() {}

  Bool setInterval(const Time &timeStart, const Time &timeEnd) override;
  UInt parameterCount()          const override {return countAParameter;}
  UInt gravityParameterCount()   const override {return gravityCount;}
  UInt rightSideCount()          const override {return rhs.size();}
  UInt arcCount()                const override {return countArc;}
  void parameterName(std::vector<ParameterName> &name) const override;

  Arc computeArc(UInt arcNo, CovarianceSstPtr covSst, CovariancePodPtr covPod1, CovariancePodPtr covPod2) override;
};

/***********************************************/

#endif /* __GROOPS__ */
