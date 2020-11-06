/***********************************************/
/**
* @file parameterNamesAcceleration.h
*
* @brief Parameter names of satellites forces.
*
* @author Torsten Mayer-Guerr
* @date 2020-05-29
*/
/***********************************************/

#ifndef __GROOPS_PARAMETERNAMESACCELERATION__
#define __GROOPS_PARAMETERNAMESACCELERATION__

// Latex documentation
static const char *docstringParameterNamesAcceleration = R"(
\subsection{Acceleration}
Parameter names of satellite acceleration \configClass{parametrization}{parametrizationAccelerationType}.
Arc related parameters are appended if an \configFile{inputfileInstrument}{instrument} is provided which
defines the arc structure.
An additional \config{object} name can be included in the parameter names.
)";

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "files/fileInstrument.h"
#include "classes/parametrizationAcceleration/parametrizationAcceleration.h"
#include "classes/parameterNames/parameterNames.h"

/***** CLASS ***********************************/

/** @brief Parameter names of satellites forces.
* @ingroup parameterNamesGroup
* @see ParameterNames */
class ParameterNamesAcceleration : public ParameterNamesBase
{
public:
  ParameterNamesAcceleration(Config &config)
  {
    try
    {
      std::string object;
      ParametrizationAccelerationPtr parameterization;
      FileName fileNameInstrument;

      readConfig(config, "object",              object,             Config::OPTIONAL, "", "object these parameters refers to, e.g. graceA, G023");
      readConfig(config, "parameterization",    parameterization,   Config::MUSTSET,  "", "");
      readConfig(config, "inputfileInstrument", fileNameInstrument, Config::OPTIONAL, "", "defines the arc structure for arc related parameters");
      if(isCreateSchema(config)) return;

      parameterization->parameterName(names);

      // arc related parameters
      InstrumentFile instrumentFile(fileNameInstrument);
      for(UInt arcNo=0; arcNo<instrumentFile.arcCount(); arcNo++)
      {
        const Arc arc = instrumentFile.readArc(arcNo);
        parameterization->setIntervalArc(arc.at(0).time, arc.back().time+medianSampling(arc.times()));

        std::vector<ParameterName> namesArc;
        parameterization->parameterNameArc(namesArc);

        const std::string str = "arc"+arcNo%"%i"s+".";
        std::for_each(namesArc.begin(), namesArc.end(), [&](auto &x) {x.type = str+x.type;});
        names.insert(names.end(), namesArc.begin(), namesArc.end());
      }

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
