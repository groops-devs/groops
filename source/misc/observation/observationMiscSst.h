/***********************************************/
/**
* @file observationMiscSst.h
*
* @brief Satellite to satellite tracking.
*
* @author Torsten Mayer-Guerr
* @date 2015-06-02
*
*/
/***********************************************/

#ifndef __GROOPS_OBSERVATIONMSICSST__
#define __GROOPS_OBSERVATIONMSICSST__

/***********************************************/

#include "base/parameterName.h"
#include "files/fileInstrument.h"
#include "classes/observation/observation.h"
#include "misc/observation/covarianceSst.h"
#include "misc/observation/covariancePod.h"

/***** TYPES ***********************************/

class ObservationMiscSst;
typedef std::shared_ptr<ObservationMiscSst> ObservationMiscSstPtr;

/***** CLASS ***********************************/

/** @brief Satellite to satellite tracking.
* @ingroup miscGroup
* @see Observation */
class ObservationMiscSst : public Observation
{
public:
  class Arc
  {
  public:
    Matrix   l, A, B;
    std::array<std::vector<Time>, 3> times; // SST, POD1, POD2
    OrbitArc pod1, pod2;
  };

  virtual ~ObservationMiscSst() {}

  virtual Arc computeArc(UInt arcNo,
                         CovarianceSstPtr covSst =nullptr,
                         CovariancePodPtr covPod1=nullptr,
                         CovariancePodPtr covPod2=nullptr) = 0;

  void observation(UInt arcNo, Matrix &l, Matrix &A, Matrix &B) override;

  /** @brief creates an derived instance of this class. */
  static ObservationMiscSstPtr create(Config &config, const std::string &name);
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class ObservationMiscSst.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a observation is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] observation Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates ObservationMiscSst */
template<> Bool readConfig(Config &config, const std::string &name, ObservationMiscSstPtr &observation, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/***********************************************/

#endif /* __GROOPS__ */
