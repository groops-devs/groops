/***********************************************/
/**
* @file observationMiscSst.cpp
*
* @brief Satellite to satellite tracking.
*
* @author Torsten Mayer-Guerr
* @date 2015-06-02
*
*/
/***********************************************/

#include "base/import.h"
#include "files/fileInstrument.h"
#include "misc/observation/covariancePod.h"
#include "misc/observation/observationMiscSstIntegral.h"
#include "misc/observation/observationMiscSstVariational.h"
#include "misc/observation/observationMiscSst.h"

/***** FUNCTIONS *******************************/

template<> Bool readConfig(Config &config, const std::string &name, ObservationMiscSstPtr &observation, Config::Appearance mustSet, const std::string &/*defaultValue*/, const std::string &/*annotation*/)
{
  try
  {
//     if(isCreateSchema(config))
//     {
//       config.xselement(name, ObservationMiscSst::typeName(), mustSet, Config::ONCE, "", annotation);
//       return FALSE;
//     }

    if(!hasName(config, name, mustSet))
      return FALSE;
    observation = ObservationMiscSst::create(config, name);
    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

ObservationMiscSstPtr ObservationMiscSst::create(Config &config, const std::string &name)
{
  try
  {
    ObservationMiscSstPtr observation;
    std::string type;
    readConfigChoice(config, name, type, Config::MUSTSET, "", "obervation equations (Sst)");
    if(readConfigChoiceElement(config, "sstIntegral",    type, "integral approach"))
      observation = ObservationMiscSstPtr(new ObservationMiscSstIntegral(config));
    if(readConfigChoiceElement(config, "sstVariational", type, "variational equations"))
      observation = ObservationMiscSstPtr(new ObservationMiscSstVariational(config));
    endChoice(config);

    return observation;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ObservationMiscSst::observation(UInt /*arcNo*/, Matrix &/*l*/, Matrix &/*A*/, Matrix &/*B*/)
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
