/***********************************************/
/**
* @file parameterNamesName.h
*
* @brief Single parameter name.
*
* @author Torsten Mayer-Guerr
* @date 2020-05-29
*/
/***********************************************/

#ifndef __GROOPS_PARAMETERNAMESNAME__
#define __GROOPS_PARAMETERNAMESNAME__

// Latex documentation
static const char *docstringParameterNamesName = R"(
\subsection{Name}
The parameter is given by explicitly by four parts:
\begin{enumerate}
\item object: Object this parameter refers to, e.g. \verb|graceA|, \verb|G023|, \verb|earth|, \ldots
\item type: Type of this parameter, e.g. \verb|accBias|, \verb|position.x|, \ldots
\item temporal: Temporal representation of this parameter, e.g. \verb|trend|, \verb|polynomial.degree1|, \ldots
\item interval: Interval/epoch this parameter represents, e.g. \verb|2017-01-01_00-00-00_2017-01-02_00-00-00|, \verb|2018-01-01_00-00-00|.
\end{enumerate}
)";

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "classes/parameterNames/parameterNames.h"

/***** CLASS ***********************************/

/** @brief Single parameter name.
* @ingroup parameterNamesGroup
* @see ParameterNames */
class ParameterNamesName : public ParameterNamesBase
{
public:
  ParameterNamesName(Config &config)
  {
    try
    {
      names.resize(1);
      readConfig(config, "object",   names.at(0).object,   Config::OPTIONAL,  "",  "object this parameter refers to, e.g. graceA, G023, earth");
      readConfig(config, "type",     names.at(0).type,     Config::OPTIONAL,  "",  "type of this parameter, e.g. accBias, position.x");
      readConfig(config, "temporal", names.at(0).temporal, Config::OPTIONAL,  "",  "temporal representation of this parameter, e.g. trend, polynomial.degree1");
      readConfig(config, "interval", names.at(0).interval, Config::OPTIONAL,  "",  "interval/epoch this parameter refers to, e.g. 2017-01-01_00-00-00_2017-01-02_00-00-00, 2008-01-01_00-00-00");
    }
    catch(std::exception &e)
    {
      GROOPS_RETHROW(e)
    }
  }
};

/***********************************************/

#endif /* __GROOPS__ */
