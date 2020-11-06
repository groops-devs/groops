/***********************************************/
/**
* @file parameterNamesGnssAntenna.h
*
* @brief Parameter names of GNSS antenna parametrization.
*
* @author Torsten Mayer-Guerr
* @date 2020-05-29
*/
/***********************************************/

#ifndef __GROOPS_PARAMETERNAMESPARAMETRIZATIONGNSSANTENNA__
#define __GROOPS_PARAMETERNAMESPARAMETRIZATIONGNSSANTENNA__

// Latex documentation
static const char *docstringParameterNamesGnssAntenna = R"(
\subsection{GnssAntenna}
Parameter names of GNSS antenna center variation \configClass{parametrization}{parametrizationGnssAntennaType}.
An additional \config{object} name (antenna name) can be included in the parameter names.
It is possible to setup the parameters for each \configClass{gnssType}{gnssType}.
)";

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "classes/parametrizationGnssAntenna/parametrizationGnssAntenna.h"
#include "classes/parameterNames/parameterNames.h"

/***** CLASS ***********************************/

/** @brief Parameter names of ParametrizationGnssAntenna.
* @ingroup parameterNamesGroup
* @see ParameterNames */
class ParameterNamesGnssAntenna : public ParameterNamesBase
{
public:
  ParameterNamesGnssAntenna(Config &config)
  {
    try
    {
      std::string                   object;
      ParametrizationGnssAntennaPtr parametrization;
      std::vector<GnssType>         types;

      readConfig(config, "object",          object,          Config::OPTIONAL, "", "antenna name");
      readConfig(config, "parametrization", parametrization, Config::MUSTSET,  "", "");
      readConfig(config, "gnssType",        types,           Config::OPTIONAL, "", "e.g. C1CG**");
      if(isCreateSchema(config)) return;

      std::vector<ParameterName> baseNames;
      parametrization->parameterName(baseNames);
      if(types.size())
        for(GnssType type : types)
        {
          std::string typeStr = "." + type.str();
          for(const auto &base : baseNames)
            names.push_back(ParameterName(object, base.type+typeStr, base.temporal, base.interval));
        }
      else
        for(const auto &base : baseNames)
          names.push_back(ParameterName(object, base.type, base.temporal, base.interval));
    }
    catch(std::exception &e)
    {
      GROOPS_RETHROW(e)
    }
  }
};

/***********************************************/

#endif /* __GROOPS__ */
