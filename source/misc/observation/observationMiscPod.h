/***********************************************/
/**
* @file observationMiscPod.h
*
* @brief Precise Orbit data.
*
* @author Torsten Mayer-Guerr
* @date 2015-06-02
*
*/
/***********************************************/

#ifndef __GROOPS_OBSERVATIONMSICPOD__
#define __GROOPS_OBSERVATIONMSICPOD__

/***********************************************/

#include "base/parameterName.h"
#include "files/fileInstrument.h"
#include "misc/observation/covariancePod.h"

/***** TYPES ***********************************/

class ObservationMiscPod;
typedef std::shared_ptr<ObservationMiscPod> ObservationMiscPodPtr;

/***** CLASS ***********************************/

/** @brief Precise Orbit data.
* @ingroup miscGroup
* @see Observation */
class ObservationMiscPod
{
public:
  class Arc
  {
  public:
    Matrix l, A, B;
    std::vector<Time> times;
    OrbitArc pod;
  };

  virtual ~ObservationMiscPod() {}

  virtual void setInterval(const Time &timeStart, const Time &timeEnd) = 0;
  virtual UInt parameterCount()          const = 0;
  virtual UInt gravityParameterCount()   const = 0;
  virtual UInt rightSideCount()          const = 0;
  virtual UInt arcCount()                const = 0;
  virtual void parameterName(std::vector<ParameterName> &name) const = 0;

  virtual Arc computeArc(UInt arcNo, CovariancePodPtr covPod=CovariancePodPtr(nullptr)) = 0;

  /** @brief creates an derived instance of this class. */
  static ObservationMiscPodPtr create(Config &config, const std::string &name);
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class ObservationMiscPod.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a observation is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] observation Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates ObservationMiscPod */
template<> Bool readConfig(Config &config, const std::string &name, ObservationMiscPodPtr &observation, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/***********************************************/

#endif /* __GROOPS__ */
