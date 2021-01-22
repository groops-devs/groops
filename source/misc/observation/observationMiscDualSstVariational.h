/***********************************************/
/**
* @file observationMiscDualSstVariational.h
*
* @brief Satellite to satellite tracking (Variational equations).
*
* @author Andreas Kvas
* @date 2019-11-22
*
*/
/***********************************************/

#ifndef __GROOPS_OBSERVATIONMSICDUALSSTVARIATIONAL__
#define __GROOPS_OBSERVATIONMSICDUALSSTVARIATIONAL__

/***********************************************/

#include "base/polynomial.h"
#include "files/fileInstrument.h"
#include "classes/observation/observation.h"
#include "classes/ephemerides/ephemerides.h"
#include "classes/parametrizationGravity/parametrizationGravity.h"
#include "classes/parametrizationAcceleration/parametrizationAcceleration.h"
#include "classes/parametrizationSatelliteTracking/parametrizationSatelliteTracking.h"
#include "misc/observation/variationalEquationFromFile.h"
#include "misc/observation/covarianceSst.h"
#include "misc/observation/covariancePod.h"

/***** TYPES ***********************************/

class ObservationMiscDualSstVariational;
typedef std::shared_ptr<ObservationMiscDualSstVariational> ObservationMiscDualSstVariationalPtr;

/***** CLASS ***********************************/

/** @brief Satellite to satellite tracking (Variational equations).
* @ingroup miscGroup
* @see Observation */
class ObservationMiscDualSstVariational : public Observation
{
public:
  class Arc
  {
  public:
    Matrix l, A, B;
    std::vector<Time>     timesSst;
    std::vector<Time>     timesPod1, timesPod2;
    OrbitArc              pod1, pod2;
    std::vector<Rotary3d> rotSat1, rotSat2;
    Matrix                pos1, pos2;
  };

  std::vector<InstrumentFilePtr>      sst1File, sst2File;
  VariationalEquationFromFile         variationalEquation1, variationalEquation2;
  InstrumentFile                      pod1File, pod2File;
  Polynomial                          polynomial;
  UInt                                countArc;
  UInt                                sstType; // 0: biased range, 1: range-rate
  Bool                                computeVelocity;
  EphemeridesPtr                      ephemerides;
  ParametrizationGravityPtr           parameterGravity;
  ParametrizationAccelerationPtr      parameterAcceleration1, parameterAcceleration2;
  ParametrizationSatelliteTrackingPtr parameterSst1, parameterSst2;

  // Indicies for design matrix A
  UInt countAParameter;
  UInt idxGravity,      gravityCount;
  UInt idxState1,       state1Count;
  UInt idxState2,       state2Count;
  UInt idxSstPara1,     idxSstPara2;

public:
  ObservationMiscDualSstVariational(Config &config);
 ~ObservationMiscDualSstVariational() {}

  Bool setInterval(const Time &timeStart, const Time &timeEnd) override;
  UInt parameterCount()          const override {return countAParameter;}
  UInt gravityParameterCount()   const override {return gravityCount;}
  UInt rightSideCount()          const override {return 1;}
  UInt arcCount()                const override {return countArc;}
  void parameterName(std::vector<ParameterName> &name) const override;

  Arc computeArc(UInt arcNo, CovarianceSstPtr covSst1=nullptr, CovarianceSstPtr covSst2=nullptr, CovarianceSstPtr covAcc=nullptr,
                 CovariancePodPtr covPod1=nullptr, CovariancePodPtr covPod2=nullptr,
                 std::vector<Rotary3d> rotSat1={}, std::vector<Rotary3d> rotSat2={});

  void observation(UInt arcNo, Matrix &l, Matrix &A, Matrix &B) override;

  static Matrix decorrelate(const_MatrixSliceRef covSst1, const_MatrixSliceRef covSst2, MatrixSliceRef covAcc, const std::list<MatrixSlice> &A);

   /** @brief creates an derived instance of this class. */
  static ObservationMiscDualSstVariationalPtr create(Config &config, const std::string &name);
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class ObservationMiscDualSst.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a observation is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] observation Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates ObservationMiscDualSst */
template<> Bool readConfig(Config &config, const std::string &name, ObservationMiscDualSstVariationalPtr &observation, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/***********************************************/

#endif /* __GROOPS__ */
