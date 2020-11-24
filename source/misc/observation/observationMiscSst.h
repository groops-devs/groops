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
#include "misc/observation/covarianceSst.h"
#include "misc/observation/covariancePod.h"

/***** TYPES ***********************************/

class ObservationMiscSst;
typedef std::shared_ptr<ObservationMiscSst> ObservationMiscSstPtr;

/***** CLASS ***********************************/

/** @brief Satellite to satellite tracking.
* @ingroup miscGroup
* @see Observation */
class ObservationMiscSst
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

  virtual ~ObservationMiscSst() {}

  virtual Bool setInterval(const Time &timeStart, const Time &timeEnd) = 0;
  virtual UInt parameterCount()          const = 0;
  virtual UInt gravityParameterCount()   const = 0;
  virtual UInt rightSideCount()          const = 0;
  virtual UInt arcCount()                const = 0;
  virtual void parameterName(std::vector<ParameterName> &name) const = 0;

  virtual Arc computeArc(UInt arcNo,
                         CovarianceSstPtr covSst =CovarianceSstPtr(nullptr),
                         CovariancePodPtr covPod1=CovariancePodPtr(nullptr),
                         CovariancePodPtr covPod2=CovariancePodPtr(nullptr),
                         const std::vector<Rotary3d> &rotSat1={},
                         const std::vector<Rotary3d> &rotSat2={}) = 0;

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
