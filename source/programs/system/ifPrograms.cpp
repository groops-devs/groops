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
Runs a list of \config{program}s if a \configClass{condition}{conditionType} is met.
Otherwise \config{elseProgram}s are executed.
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
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(IfPrograms, PARALLEL, "Runs programs if condition is met.", System)
GROOPS_RENAMED_PROGRAM(IfProgramme, IfPrograms, date2time(2020, 6, 3))

/***********************************************/

void IfPrograms::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    ConditionPtr  conditionPtr;
    ProgramConfig programs, elsePrograms;

    renameDeprecatedConfig(config, "programme", "program", date2time(2020, 6, 3));

    readConfig(config, "condition",   conditionPtr, Config::MUSTSET,  "", "");
    readConfig(config, "program",     programs,     Config::OPTIONAL, "", "executed if condition evaluates to true");
    readConfig(config, "elseProgram", elsePrograms, Config::OPTIONAL, "", "executed if condition evaluates to false");
    if(isCreateSchema(config)) return;

    VariableList varList;
    if(conditionPtr->condition(varList))
    {
      logInfo<<"  condition is true."<<Log::endl;
      programs.run(varList, comm);
    }
    else
    {
      logInfo<<"  condition is false."<<Log::endl;
      elsePrograms.run(varList, comm);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
