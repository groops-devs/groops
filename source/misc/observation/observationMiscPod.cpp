/***********************************************/
/**
* @file observationMiscPod.cpp
*
* @brief Precise Orbit data.
*
* @author Torsten Mayer-Guerr
* @date 2015-06-02
*
*/
/***********************************************/

#include "base/import.h"
#include "files/fileInstrument.h"
#include "misc/observation/covariancePod.h"
#include "misc/observation/observationMiscPodIntegral.h"
#include "misc/observation/observationMiscPodVariational.h"
#include "misc/observation/observationMiscPod.h"

/***** FUNCTIONS *******************************/

template<> Bool readConfig(Config &config, const std::string &name, ObservationMiscPodPtr &observation, Config::Appearance mustSet, const std::string &/*defaultValue*/, const std::string &/*annotation*/)
{
  try
  {
//     if(isCreateSchema(config))
//     {
//       config.xselement(name, ObservationMiscPod::typeName(), mustSet, Config::ONCE, "", annotation);
//       return FALSE;
//     }

    if(!hasName(config, name, mustSet))
      return FALSE;
    observation = ObservationMiscPod::create(config, name);
    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

ObservationMiscPodPtr ObservationMiscPod::create(Config &config, const std::string &name)
{
  try
  {
    ObservationMiscPodPtr observation;
    std::string type;
    readConfigChoice(config, name, type, Config::MUSTSET, "", "obervation equations (POD)");
    if(readConfigChoiceElement(config, "podIntegral"        , type, "precise orbit data (integral approach)"))
      observation = ObservationMiscPodPtr(new ObservationMiscPodIntegral(config));
    if(readConfigChoiceElement(config, "podVariational"     , type, "precise orbit data (variational equations)"))
      observation = ObservationMiscPodPtr(new ObservationMiscPodVariational(config));
    endChoice(config);

    return observation;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ObservationMiscPod::observation(UInt /*arcNo*/, Matrix &/*l*/, Matrix &/*A*/, Matrix &/*B*/)
{
  try
  {
    throw(Exception("Must not be called"));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
