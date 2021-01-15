/***********************************************/
/**
* @file groupPrograms.cpp
*
* @brief Runs programs in group and can catch errors.
*
* @author Sebastian Strasser
* @author Torsten Mayer-Guerr
* @date 2017-09-01
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Runs \config{program}s in a group, which can be used to structure a config file.
If \config{catchErrors} is enabled and an error occurs, the remaining \config{program}s
are skipped and execution continues with \config{errorProgram}s, in case any are defined.
Otherwise an exception is thrown.
)";

/***********************************************/

#include "programs/program.h"

/***** CLASS ***********************************/

/** @brief  Runs programs in group and can catch errors.
* @ingroup programsGroup */
class GroupPrograms
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GroupPrograms, PARALLEL, "Runs programs in group and can catch errors.", System)
GROOPS_RENAMED_PROGRAM(GroupProgramme, GroupPrograms, date2time(2020, 6, 3))

/***********************************************/

void GroupPrograms::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    Bool          continueAfterError = FALSE;
    ProgramConfig programs, errorPrograms;

    renameDeprecatedConfig(config, "programme", "program", date2time(2020, 6, 3));

    readConfig(config, "program", programs, Config::OPTIONAL, "", "");
    if(readConfigSequence(config, "catchErrors", Config::OPTIONAL, "", ""))
    {
      continueAfterError = TRUE;
      readConfig(config, "errorProgram", errorPrograms, Config::OPTIONAL, "", "executed if an error occured");
      endSequence(config);
    }
    if(isCreateSchema(config)) return;

    try
    {
      Parallel::broadCastExceptions(comm, [&](Parallel::CommunicatorPtr comm)
      {
        auto varListTmp = config.getVarList();
        programs.run(varListTmp, comm);
      });
    }
    catch(std::exception &e)
    {
      if(!continueAfterError)
        throw;
      if(Parallel::isMaster(comm))
        logError<<e.what()<<"  continue with error programs..."<<Log::endl;
      auto varList = config.getVarList();
      errorPrograms.run(varList, comm);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
