/***********************************************/
/**
* @file parameterNamesTemporal.h
*
* @brief Parameter names from temporal parametrization.
*
* @author Torsten Mayer-Guerr
* @date 2020-05-29
*/
/***********************************************/

#ifndef __GROOPS_PARAMETERNAMESTEMPORAL__
#define __GROOPS_PARAMETERNAMESTEMPORAL__

// Latex documentation
static const char *docstringParameterNamesTemporal = R"(
\subsection{Temporal}
Parameter names from temporal parametrization.
It is possible to setup the temporal parameters for each \configClass{parameterNameBase}{parameterNamesType}.
)";

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "classes/parameterNames/parameterNames.h"

/***** CLASS ***********************************/

/** @brief Parameter names from temporal parametrization.
* @ingroup parameterNamesGroup
* @see ParameterNames */
class ParameterNamesTemporal : public ParameterNamesBase
{
public:
  ParameterNamesTemporal(Config &config)
  {
    try
    {
      ParameterNamesPtr          parameterNames;
      ParametrizationTemporalPtr temporal;

      readConfig(config, "parameterNameBase",       parameterNames, Config::OPTIONAL, "", "");
      readConfig(config, "parametrizationTemporal", temporal,       Config::MUSTSET,  "", "");
      if(isCreateSchema(config)) return;

      if(parameterNames && parameterNames->parameterNames().size())
        temporal->parameterName(parameterNames->parameterNames(), names);
      else
        temporal->parameterName(names);
    }
    catch(std::exception &e)
    {
      GROOPS_RETHROW(e)
    }
  }
};

/***********************************************/

#endif /* __GROOPS__ */
