/***********************************************/
/**
* @file parameterNamesRename.h
*
* @brief Single parameter name.
*
* @author Torsten Mayer-Guerr
* @date 2020-05-29
*/
/***********************************************/

#ifndef __GROOPS_PARAMETERNAMESRENAME__
#define __GROOPS_PARAMETERNAMESRENAME__

// Latex documentation
static const char *docstringParameterNamesRename = R"(
\subsection{Rename}
Replaces parts of \configClass{parameterName}{parameterNamesType}s.
The star "\verb|*|" left this part untouched.
)";

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "classes/parameterNames/parameterNames.h"

/***** CLASS ***********************************/

/** @brief Single parameter name.
* @ingroup parameterNamesGroup
* @see ParameterNames */
class ParameterNamesRename : public ParameterNamesBase
{
public:
  ParameterNamesRename(Config &config)
  {
    try
    {
      ParameterNamesPtr parameterNames;
      std::string       object, type, temporal, interval;

      readConfig(config, "parameterName", parameterNames, Config::MUSTSET,  "",  "");
      readConfig(config, "object",        object,         Config::OPTIONAL, "*", "*: left this part untouched, object");
      readConfig(config, "type",          type,           Config::OPTIONAL, "*", "*: left this part untouched, type");
      readConfig(config, "temporal",      temporal,       Config::OPTIONAL, "*", "*: left this part untouched, temporal representation");
      readConfig(config, "interval",      interval,       Config::OPTIONAL, "*", "*: left this part untouched, interval/epoch");
      if(isCreateSchema(config)) return;

      names = parameterNames->parameterNames();
      for(auto &name : names)
      {
        if(object   != "*") name.object   = object;
        if(type     != "*") name.type     = type;
        if(temporal != "*") name.temporal = temporal;
        if(interval != "*") name.interval = interval;
      }
    }
    catch(std::exception &e)
    {
      GROOPS_RETHROW(e)
    }
  }
};

/***********************************************/

#endif /* __GROOPS__ */
