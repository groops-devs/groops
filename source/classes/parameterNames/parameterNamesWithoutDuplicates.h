/***********************************************/
/**
* @file parameterNamesWithoutDuplicates.h
*
* @brief Removes all duplicate names.
*
* @author Torsten Mayer-Guerr
* @date 2020-05-29
*/
/***********************************************/

#ifndef __GROOPS_PARAMETERNAMESWITHOUTDUPLICATES__
#define __GROOPS_PARAMETERNAMESWITHOUTDUPLICATES__

// Latex documentation
static const char *docstringParameterNamesWithoutDuplicates = R"(
\subsection{WithoutDuplicates}
Removes all duplicate names (keep first) from \configClass{parameterName}{parameterNamesType}.
)";

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "classes/parameterNames/parameterNames.h"

/***** CLASS ***********************************/

/** @brief Removes all duplicate names.
* @ingroup parameterNamesGroup
* @see ParameterNames */
class ParameterNamesWithoutDuplicates : public ParameterNamesBase
{
public:
  ParameterNamesWithoutDuplicates(Config &config)
  {
    try
    {
      ParameterNamesPtr parameterNames;

      readConfig(config, "parameterName", parameterNames, Config::MUSTSET,  "", "");
      if(isCreateSchema(config)) return;

      std::set<ParameterName> sortedNames;
      for(const auto &parameterName : parameterNames->parameterNames())
        if(sortedNames.find(parameterName) == sortedNames.end())
        {
          sortedNames.insert(parameterName);
          names.push_back(parameterName);
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
