/***********************************************/
/**
* @file parameterNamesCreate.cpp
*
* @brief Generate a parameter names file.
*
* @author Torsten Mayer-Guerr
* @date 2019-09-28
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Generate a \configFile{outputfileParameterNames}{parameterName} by \configClass{parameterName}{parameterNamesType}.
This file can be used in \program{NormalsCreate} or in the class \configClass{parameterSelector}{parameterSelectorType}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileParameterName.h"
#include "classes/parameterNames/parameterNames.h"

/***** CLASS ***********************************/

/** @brief Generate a parameter names file.
* @ingroup programsGroup */
class ParameterNamesCreate
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(ParameterNamesCreate, SINGLEPROCESS, "Generate a parameter names file.", NormalEquation)

/***********************************************/

void ParameterNamesCreate::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName          fileNameParameterNames;
    ParameterNamesPtr parameterNames;

    readConfig(config, "outputfileParameterNames", fileNameParameterNames, Config::MUSTSET, "", "output parameter names file");
    readConfig(config, "parameterName",            parameterNames,         Config::MUSTSET, "", "");
    if(isCreateSchema(config)) return;

    logStatus<<"writing parameter names file <"<<fileNameParameterNames<<">"<<Log::endl;
    writeFileParameterName(fileNameParameterNames, parameterNames->parameterNames());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
