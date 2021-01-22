/***********************************************/
/**
* @file observationMiscSstVariational.h
*
* @brief Satellite to satellite tracking (Variational equations).
*
* @author Torsten Mayer-Guerr
* @date 2015-06-02
*
*/
/***********************************************/

#ifndef __GROOPS_OBSERVATIONMSICSSTVARIATIONAL__
#define __GROOPS_OBSERVATIONMSICSSTVARIATIONAL__

/***********************************************/

#include "base/polynomial.h"
#include "files/fileInstrument.h"
#include "classes/ephemerides/ephemerides.h"
#include "classes/parametrizationGravity/parametrizationGravity.h"
#include "classes/parametrizationAcceleration/parametrizationAcceleration.h"
#include "classes/parametrizationSatelliteTracking/parametrizationSatelliteTracking.h"
#include "misc/observation/variationalEquationFromFile.h"
#include "misc/observation/observationMiscSst.h"

/***** TYPES ***********************************/

class ObservationMiscSstVariational;
typedef std::shared_ptr<ObservationMiscSstVariational> ObservationMiscSstVariationalPtr;

/***** CLASS ***********************************/

/** @brief Satellite to satellite tracking (Variational equations).
* @ingroup miscGroup
* @see Observation */
class ObservationMiscSstVariational : public ObservationMiscSst
{
public:
  std::vector<InstrumentFilePtr>      sstFile;
  VariationalEquationFromFile         variationalEquation1, variationalEquation2;
  InstrumentFile                      pod1File, pod2File;
  Polynomial                          polynomial;
  UInt                                countArc;
  UInt                                sstType; // 0: biased range, 1: range-rate
  Bool                                computeVelocity;
  EphemeridesPtr                      ephemerides;
  ParametrizationGravityPtr           parameterGravity;
  ParametrizationAccelerationPtr      parameterAcceleration1, parameterAcceleration2;
  ParametrizationSatelliteTrackingPtr parameterSst;

  // Indicies for design matrix A
  UInt countAParameter;
  UInt idxGravity,      gravityCount;
  UInt idxState1,       state1Count;
  UInt idxState2,       state2Count;
  UInt idxSstPara;

public:
  ObservationMiscSstVariational(Config &config);
 ~ObservationMiscSstVariational() {}

  Bool setInterval(const Time &timeStart, const Time &timeEnd) override;
  UInt parameterCount()          const override {return countAParameter;}
  UInt gravityParameterCount()   const override {return gravityCount;}
  UInt rightSideCount()          const override {return 1;}
  UInt arcCount()                const override {return countArc;}
  void parameterName(std::vector<ParameterName> &name) const override;

  Arc computeArc(UInt arcNo, CovarianceSstPtr covSst,
                 CovariancePodPtr covPod1, CovariancePodPtr covPod2,
                 const std::vector<Rotary3d> &rotSat1={}, const std::vector<Rotary3d> &rotSat2={}) override;
};

/***********************************************/

#endif /* __GROOPS__ */
