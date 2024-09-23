/***********************************************/
/**
* @file parameterNames.cpp
*
* @brief Generate parameter names.
*
* @author Torsten Mayer-Guerr
* @date 2020-05-29
*/
/***********************************************/

#define DOCSTRING_ParameterNames

#include "base/import.h"
#include "config/configRegister.h"
#include "classes/parameterNames/parameterNamesName.h"
#include "classes/parameterNames/parameterNamesFile.h"
#include "classes/parameterNames/parameterNamesGravity.h"
#include "classes/parameterNames/parameterNamesAcceleration.h"
#include "classes/parameterNames/parameterNamesSatelliteTracking.h"
#include "classes/parameterNames/parameterNamesTemporal.h"
#include "classes/parameterNames/parameterNamesGnssAntenna.h"
#include "classes/parameterNames/parameterNamesObservation.h"
#include "classes/parameterNames/parameterNamesRename.h"
#include "classes/parameterNames/parameterNamesSelection.h"
#include "classes/parameterNames/parameterNamesWithoutDuplicates.h"
#include "classes/parameterNames/parameterNames.h"

/***********************************************/

GROOPS_REGISTER_CLASS(ParameterNames, "parameterNamesType",
                      ParameterNamesName,
                      ParameterNamesFile,
                      ParameterNamesGravity,
                      ParameterNamesAcceleration,
                      ParameterNamesSatelliteTracking,
                      ParameterNamesTemporal,
                      ParameterNamesGnssAntenna,
                      ParameterNamesObservation,
                      ParameterNamesRename,
                      ParameterNamesSelection,
                      ParameterNamesWithoutDuplicates)

GROOPS_READCONFIG_UNBOUNDED_CLASS(ParameterNames, "parameterNamesType")

/***********************************************/

ParameterNames::ParameterNames(Config &config, const std::string &name)
{
  try
  {
    std::unique_ptr<ParameterNamesBase> ptr;

    std::string choice;
    while(readConfigChoice(config, name, choice, Config::OPTIONAL, "", "generate a list of parameter names"))
    {
      if(readConfigChoiceElement(config, "name",                             choice, ""))
        ptr = std::unique_ptr<ParameterNamesBase>(new ParameterNamesName(config));
      if(readConfigChoiceElement(config, "file",                             choice, ""))
        ptr = std::unique_ptr<ParameterNamesBase>(new ParameterNamesFile(config));
      if(readConfigChoiceElement(config, "parametrizationGravity",           choice, "parametrization of gravity field"))
        ptr = std::unique_ptr<ParameterNamesBase>(new ParameterNamesGravity(config));
      if(readConfigChoiceElement(config, "parametrizationAcceleration",      choice, "parametrization of orbit forces"))
        ptr = std::unique_ptr<ParameterNamesBase>(new ParameterNamesAcceleration(config));
      if(readConfigChoiceElement(config, "parametrizationSatelliteTracking", choice, "satellite tracking parametrization"))
        ptr = std::unique_ptr<ParameterNamesBase>(new ParameterNamesSatelliteTracking(config));
      if(readConfigChoiceElement(config, "parametrizationTemporal",          choice, "temporal parametrization"))
        ptr = std::unique_ptr<ParameterNamesBase>(new ParameterNamesTemporal(config));
      if(readConfigChoiceElement(config, "parametrizationGnssAntenna",       choice, "parametrization of GNSS antenna center variations"))
        ptr = std::unique_ptr<ParameterNamesBase>(new ParameterNamesGnssAntenna(config));
      if(readConfigChoiceElement(config, "observation",                      choice, "parameters of oberservation equations"))
        ptr = std::unique_ptr<ParameterNamesBase>(new ParameterNamesObservation(config));
      if(readConfigChoiceElement(config, "rename",                           choice, "rename parts"))
        ptr = std::unique_ptr<ParameterNamesBase>(new ParameterNamesRename(config));
      if(readConfigChoiceElement(config, "selection",                        choice, ""))
        ptr = std::unique_ptr<ParameterNamesBase>(new ParameterNamesSelection(config));
      if(readConfigChoiceElement(config, "withoutDuplicates",                choice, ""))
        ptr = std::unique_ptr<ParameterNamesBase>(new ParameterNamesWithoutDuplicates(config));
      endChoice(config);
      if(isCreateSchema(config))
        return;

      names.insert(names.end(), ptr->names.begin(), ptr->names.end());
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
