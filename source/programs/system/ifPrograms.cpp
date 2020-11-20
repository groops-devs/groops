/***********************************************/
/**
* @file ifPrograms.cpp
*
* @brief Runs programs if condition is met.
*
* @author Andreas Kvas
* @date 2017-01-04
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Runs a list of programs if a \configClass{condition}{conditionType} is met.
)";

/***********************************************/

#include "programs/program.h"
#include "classes/condition/condition.h"

/***** CLASS ***********************************/

/** @brief  Runs programs if condition is met.
* @ingroup programsGroup */
class IfPrograms
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(IfPrograms, PARALLEL, "Runs programs if condition is met.", System)
GROOPS_RENAMED_PROGRAM(IfProgramme, IfPrograms, date2time(2020, 6, 3))

/***********************************************/

void IfPrograms::run(Config &config)
{
  try
  {
    ConditionPtr  conditionPtr;
    ProgramConfig programs;

    renameDeprecatedConfig(config, "programme", "program", date2time(2020, 6, 3));

    readConfig(config, "condition", conditionPtr, Config::MUSTSET,  "", "");
    readConfig(config, "program",   programs,     Config::OPTIONAL, "", "");
    if(isCreateSchema(config)) return;

    auto varList = config.getVarList();
    if(conditionPtr->condition(varList))
    {
      logInfo<<"  condition is true."<<Log::endl;
      programs.run(varList);
    }
    else
    {
      logInfo<<"  condition is false."<<Log::endl;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
