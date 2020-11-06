/***********************************************/
/**
* @file parameterNamesGravity.h
*
* @brief Parameter names from gravity field parametrization.
*
* @author Torsten Mayer-Guerr
* @date 2020-05-29
*/
/***********************************************/

#ifndef __GROOPS_PARAMETERNAMESGRAVITY__
#define __GROOPS_PARAMETERNAMESGRAVITY__

// Latex documentation
static const char *docstringParameterNamesGravity = R"(
\subsection{Gravity}
Parameter names of gravity \configClass{parametrization}{parametrizationGravityType}.
An additional \config{object} name can be included in the parameter names.
)";

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "classes/parametrizationGravity/parametrizationGravity.h"
#include "classes/parameterNames/parameterNames.h"

/***** CLASS ***********************************/

/** @brief Parameter names from gravity field parametrization.
* @ingroup parameterNamesGroup
* @see ParameterNames */
class ParameterNamesGravity : public ParameterNamesBase
{
public:
  ParameterNamesGravity(Config &config)
  {
    try
    {
      std::string object;
      ParametrizationGravityPtr parametrization;

      readConfig(config, "object",          object,          Config::OPTIONAL, "", "object these parameters refers to, e.g. earth");
      readConfig(config, "parametrization", parametrization, Config::MUSTSET,  "", "");
      if(isCreateSchema(config)) return;

      parametrization->parameterName(names);
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
