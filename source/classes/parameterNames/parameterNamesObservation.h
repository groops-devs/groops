/***********************************************/
/**
* @file parameterNamesObservation.h
*
* @brief Parameter names from observation equations.
*
* @author Torsten Mayer-Guerr
* @date 2020-05-29
*/
/***********************************************/

#ifndef __GROOPS_PARAMETERNAMESOBSERVATION__
#define __GROOPS_PARAMETERNAMESOBSERVATION__

// Latex documentation
static const char *docstringParameterNamesObservation = R"(
\subsection{Observation}
Parameter names used in \configClass{observation equations}{observationType}.
)";

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "classes/observation/observation.h"
#include "classes/parameterNames/parameterNames.h"

/***** CLASS ***********************************/

/** @brief Parameter names from observation equations.
* @ingroup parameterNamesGroup
* @see ParameterNames */
class ParameterNamesObservation : public ParameterNamesBase
{
public:
  ParameterNamesObservation(Config &config)
  {
    try
    {
      ObservationPtr observation;

      readConfig(config, "observation", observation, Config::MUSTSET,  "",  "");
      if(isCreateSchema(config)) return;

      observation->parameterName(names);
    }
    catch(std::exception &e)
    {
      GROOPS_RETHROW(e)
    }
  }
};

/***********************************************/

#endif /* __GROOPS__ */
