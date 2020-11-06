/***********************************************/
/**
* @file parameterSelectorWildcard.h
*
* @brief Parameter index vector from name.
*
* @author Sebastian Strasser
* @date 2018-05-08
*
*/
/***********************************************/

#ifndef __GROOPS_PARAMETERSELECTORNAME__
#define __GROOPS_PARAMETERSELECTORNAME__

// Latex documentation
#ifdef DOCSTRING_ParameterSelector
static const char *docstringParameterSelectorWildcard = R"(
\subsection{Wildcard}\label{parameterSelectorType:wildcard}
Parameter index vector from name. Name matching supports wildcards * for any number of characters and ? for exactly one character.
Does not add zero/empty parameters if there are no matches.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "base/string.h"
#include "classes/parameterSelector/parameterSelector.h"

/***** CLASS ***********************************/

/** @brief Parameter index vector from name.
* @ingroup parameterSelectorGroup
* @see ParameterSelector */
class ParameterSelectorWildcard : public ParameterSelectorBase
{
  std::string object, type, temporal, interval;

public:
  ParameterSelectorWildcard(Config &config);
  std::vector<UInt> indexVector(const std::vector<ParameterName> &parameterWildcards, VariableList varList);
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline ParameterSelectorWildcard::ParameterSelectorWildcard(Config &config)
{
  try
  {
    readConfig(config, "object",   object,   Config::OPTIONAL, "*", "object this parameter refers to, e.g. graceA, G023, earth (wildcards: * and ?)");
    readConfig(config, "type",     type,     Config::OPTIONAL, "*", "type of this parameter, e.g. accBias, position.x (wildcards: * and ?)");
    readConfig(config, "temporal", temporal, Config::OPTIONAL, "*", "temporal representation of this parameter, e.g. trend, polynomial.degree1 (wildcards: * and ?)");
    readConfig(config, "interval", interval, Config::OPTIONAL, "*", "interval/epoch this parameter refers to, e.g. 2017-01-01_00-00-00_2017-01-02_00-00-00, 2008-01-01_00-00-00 (wildcards: * and ?)");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline std::vector<UInt> ParameterSelectorWildcard::indexVector(const std::vector<ParameterName> &parameterWildcards, VariableList /*varList*/)
{
  try
  {
    const std::regex pattern = String::wildcard2regex(ParameterName(object, type, temporal, interval).str());

    std::vector<UInt> vector;
    for(UInt i=0; i<parameterWildcards.size(); i++)
      if(std::regex_match(parameterWildcards.at(i).str(), pattern))
        vector.push_back(i);

    return vector;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
