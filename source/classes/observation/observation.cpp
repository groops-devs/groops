/***********************************************/
/**
* @file observation.cpp
*
* @brief Observation equations.
* Set up linearized observation equations (design matrix)
* to connect unknown parameters with observations.
* It is used in an least squares adjustement.
*
* @author Torsten Mayer-Guerr
* @date 2001-11-12
*
*/
/***********************************************/

#define DOCSTRING_Observation

#include "base/import.h"
#include "config/configRegister.h"
#include "parallel/parallel.h"
#include "observationPodIntegral.h"
#include "observationPodAcceleration.h"
#include "observationPodEnergy.h"
#include "observationPodVariational.h"
#include "observationSstIntegral.h"
#include "observationSstVariational.h"
#include "observationDualSstVariational.h"
#include "observationGradiometer.h"
#include "observationTerrestrial.h"
#include "observationDeflections.h"
#include "observationStationLoading.h"
#include "observation.h"

/***********************************************/

GROOPS_REGISTER_CLASS(Observation, "observationType",
                      ObservationPodVariational,
                      ObservationPodIntegral,
                      ObservationPodAcceleration,
                      ObservationPodEnergy,
                      ObservationSstVariational,
                      ObservationSstIntegral,
                      ObservationDualSstVariational,
                      ObservationGradiometer,
                      ObservationTerrestrial,
                      ObservationDeflections,
                      ObservationStationLoading)

GROOPS_READCONFIG_CLASS(Observation, "observationType")

/***********************************************/

ObservationPtr Observation::create(Config &config, const std::string &name)
{
  try
  {
    ObservationPtr ptr;
    std::string    type;

    readConfigChoice(config, name, type, Config::MUSTSET, "", "obervation equations for least squares adjustment");

    renameDeprecatedChoice(config, type, "satelliteTracking", "sstIntegral", date2time(2020, 6, 13));

    if(readConfigChoiceElement(config, "podVariational",      type, "precise orbit data (variational equations)"))
      ptr = std::make_shared<ObservationPodVariational>(config);
    if(readConfigChoiceElement(config, "podIntegral",         type, "precise orbit data (integral approach)"))
      ptr = std::make_shared<ObservationPodIntegral>(config);
    if(readConfigChoiceElement(config, "podAcceleration",     type, "precise orbit data (acceleration approach)"))
      ptr = std::make_shared<ObservationPodAcceleration>(config);
    if(readConfigChoiceElement(config, "podEnergy",           type, "precise orbit data (energy approach)"))
      ptr = std::make_shared<ObservationPodEnergy>(config);
    if(readConfigChoiceElement(config, "sstVariational",      type, "sst data and pod"))
      ptr = std::make_shared<ObservationSstVariational>(config);
    if(readConfigChoiceElement(config, "sstIntegral",         type, "sst data and pod"))
      ptr = std::make_shared<ObservationSstIntegral>(config);
    if(readConfigChoiceElement(config, "dualSstVariational",  type, "dual sst data and pod"))
      ptr = std::make_shared<ObservationDualSstVariational>(config);
    if(readConfigChoiceElement(config, "gradiometer",         type, "GOCE gradiometry"))
      ptr = std::make_shared<ObservationGradiometer>(config);
    if(readConfigChoiceElement(config, "terrestrial",         type, "e.g. gravity anomalies, GPS/levelling"))
      ptr = std::make_shared<ObservationTerrestrial>(config);
    if(readConfigChoiceElement(config, "deflections",         type, "e.g. Deflections of the vertical"))
      ptr = std::make_shared<ObservationDeflections>(config);
    if(readConfigChoiceElement(config, "stationLoading",      type, "Loading from station observations."))
      ptr = std::make_shared<ObservationStationLoading>(config);
    endChoice(config);

    return ptr;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
