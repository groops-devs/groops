/***********************************************/
/**
* @file parameterNamesFile.h
*
* @brief Parameter names from file.
*
* @author Torsten Mayer-Guerr
* @date 2020-05-29
*/
/***********************************************/

#ifndef __GROOPS_PARAMETERNAMESFILE__
#define __GROOPS_PARAMETERNAMESFILE__

// Latex documentation
static const char *docstringParameterNamesFile = R"(
\subsection{File}
Read parameter names from \file{file}{parameterName}.
)";

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "files/fileParameterName.h"
#include "classes/parameterNames/parameterNames.h"

/***** CLASS ***********************************/

/** @brief Parameter names from file.
* @ingroup parameterNamesGroup
* @see ParameterNames */
class ParameterNamesFile : public ParameterNamesBase
{
public:
  ParameterNamesFile(Config &config)
  {
    try
    {
      FileName fileName;

      readConfig(config, "inputfileParameterNames", fileName, Config::MUSTSET, "", "file with parameter names");
      if(isCreateSchema(config)) return;

      readFileParameterName(fileName, names);
    }
    catch(std::exception &e)
    {
      GROOPS_RETHROW(e)
    }
  }
};

/***********************************************/

#endif /* __GROOPS__ */
