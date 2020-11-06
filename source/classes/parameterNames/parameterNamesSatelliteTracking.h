/***********************************************/
/**
* @file parameterNamesSatelliteTracking.h
*
* @brief Parameter names of satellite tracking parametrization.
*
* @author Torsten Mayer-Guerr
* @date 2020-05-29
*/
/***********************************************/

#ifndef __GROOPS_PARAMETERNAMESSATELLITETRACKING__
#define __GROOPS_PARAMETERNAMESSATELLITETRACKING__

// Latex documentation
static const char *docstringParameterNamesSatelliteTracking = R"(
\subsection{SatelliteTracking}
Parameter names of satellite tracking \configClass{parametrization}{parametrizationSatelliteTrackingType}.
An additional \config{object} name can be included in the parameter names.
)";

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "classes/parametrizationSatelliteTracking/parametrizationSatelliteTracking.h"
#include "classes/parameterNames/parameterNames.h"

/***** CLASS ***********************************/

/** @brief Parameter names of satellite tracking parametrization.
* @ingroup parameterNamesGroup
* @see ParameterNames */
class ParameterNamesSatelliteTracking : public ParameterNamesBase
{
public:
  ParameterNamesSatelliteTracking(Config &config)
  {
    try
    {
      std::string object;
      ParametrizationSatelliteTrackingPtr parameterization;

      readConfig(config, "object",           object,           Config::OPTIONAL, "", "object these parameters refers to, e.g. grace1.grace2");
      readConfig(config, "parameterization", parameterization, Config::MUSTSET,  "", "");
      if(isCreateSchema(config)) return;

      parameterization->parameterName(names);
      std::for_each(names.begin(), names.end(), [&](auto &x) {x.object = object;});
    }
    catch(std::exception &e)
    {
      GROOPS_RETHROW(e)
    }
  }
};

/***********************************************/

#endif /* __GROOPS__ */
