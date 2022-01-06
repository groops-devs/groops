/***********************************************/
/**
* @file observationMiscPodVariational.h
*
* @brief Precise Orbit data (variational equations).
*
* @author Torsten Mayer-Guerr
* @date 2015-06-02
*
*/
/***********************************************/

#ifndef __GROOPS_OBSERVATIONMSICPODVARIATIONAL__
#define __GROOPS_OBSERVATIONMSICPODVARIATIONAL__

/***********************************************/

#include "files/fileInstrument.h"
#include "classes/ephemerides/ephemerides.h"
#include "classes/parametrizationGravity/parametrizationGravity.h"
#include "classes/parametrizationAcceleration/parametrizationAcceleration.h"
#include "misc/observation/variationalEquationFromFile.h"
#include "misc/observation/observationMiscPod.h"

/***** TYPES ***********************************/

class ObservationMiscPodVariational;
typedef std::shared_ptr<ObservationMiscPodVariational> ObservationMiscPodVariationalPtr;

/***** CLASS ***********************************/

/** @brief Precise Orbit data (variational equations).
* @ingroup miscGroup
* @see Observation */
class ObservationMiscPodVariational : public ObservationMiscPod
{
  InstrumentFile                 podFile;
  VariationalEquationFromFile    variationalEquation;
  UInt                           interpolationDegree;
  UInt                           countArc;
  EphemeridesPtr                 ephemerides;
  ParametrizationGravityPtr      parameterGravity;
  ParametrizationAccelerationPtr parameterAcceleration;
  Bool                           accelerateComputation;

  // Indicies for design matrix A
  UInt countAParameter;
  UInt idxGravity, gravityCount;
  UInt idxState;

public:
  ObservationMiscPodVariational(Config &config);
 ~ObservationMiscPodVariational() {}

  Bool setInterval(const Time &timeStart, const Time &timeEnd) override;
  UInt parameterCount()          const override {return countAParameter;}
  UInt gravityParameterCount()   const override {return gravityCount;}
  UInt rightSideCount()          const override {return 1;}
  UInt arcCount()                const override {return countArc;}
  void parameterName(std::vector<ParameterName> &name) const override;

  Arc computeArc(UInt arcNo, CovariancePodPtr cov) override;
};

/***********************************************/

#endif /* __GROOPS__ */
