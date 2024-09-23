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

The \config{silently} option disables the screen ouput of the \config{program}s.
With \config{outputfileLog} a log file is written for this group additional to a global log file.
This might be helpful within \program{LoopPrograms} with parallel iterations.
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
  Log::GroupPtr groupPtr;
  try {Parallel::broadCastExceptions(comm, [&](Parallel::CommunicatorPtr comm)
  {
    Bool          silently;
    ProgramConfig programs, errorPrograms;
    Bool          continueAfterError = FALSE;
    FileName      fileNameLog;

    renameDeprecatedConfig(config, "programme", "program", date2time(2020, 6, 3));

    readConfig(config, "outputfileLog", fileNameLog, Config::OPTIONAL, "",  "additional log file");
    readConfig(config, "silently",      silently,    Config::DEFAULT,  "0", "without showing the output.");
    readConfig(config, "program",       programs,    Config::OPTIONAL, "",  "");
    if(readConfigSequence(config, "catchErrors", Config::OPTIONAL, "", ""))
    {
      continueAfterError = TRUE;
      readConfig(config, "errorProgram", errorPrograms, Config::OPTIONAL, "", "executed if an error occured");
      endSequence(config);
    }
    if(isCreateSchema(config)) return;

    // -------------------------------

    Parallel::barrier(comm);
    groupPtr = Log::group(Parallel::isMaster(comm), silently);
    Log::setLogFile(fileNameLog);
    Log::currentLogFileOnly(TRUE);
    if(Parallel::size(comm) > 1)
      logStatus<<"=== Starting GROOPS subgroup with "<<Parallel::size(comm)<<" processes ==="<<Log::endl;
    else
      logStatus<<"=== Starting GROOPS subgroup ==="<<Log::endl;
    Log::currentLogFileOnly(FALSE);
    Parallel::barrier(comm);

    try
    {
      Parallel::broadCastExceptions(comm, [&](Parallel::CommunicatorPtr comm)
      {
        VariableList varList;
        programs.run(varList, comm);
      });
    }
    catch(std::exception &e)
    {
      if(!continueAfterError || Parallel::isExternal(e))
        throw;
      if(Parallel::isMaster(comm))
        logError<<e.what()<<"  continue with error programs..."<<Log::endl;
      VariableList varList;
      errorPrograms.run(varList, comm);
    }

    Parallel::barrier(comm);
    Log::currentLogFileOnly(TRUE);
    logStatus<<"=== Finished GROOPS subgroup ==="<<Log::endl;
  });}
  catch(std::exception &e)
  {
    if(groupPtr)
    {
      Log::currentLogFileOnly(TRUE);
      if(Parallel::isMaster(comm))
      {
        logError<<"****** Error ******"<<Log::endl;
        logError<<e.what()<<Log::endl;
        logStatus<<"=== Finished GROOPS subgroup with error ==="<<Log::endl;
      }
    }
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
