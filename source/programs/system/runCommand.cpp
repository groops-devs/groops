/***********************************************/
/**
* @file runCommand.cpp
*
* @brief Execute system commands.
*
* @author Matthias Ellmer
* @author Torsten Mayer-Guerr
* @date 2016-07-13
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Execute system \config{command}s. If \config{executeParallel} is set and
multiple \config{command}s are given they are executed in parallel at
distributed nodes, otherwise they are executed consecutively at master node only.
)";

/***********************************************/

#include "programs/program.h"
#include "inputOutput/system.h"

/***** CLASS ***********************************/

/** @brief Execute system commands.
* @ingroup programsGroup */
class RunCommand
{
 public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(RunCommand, PARALLEL, "Execute system commands", System)

/***********************************************/

void RunCommand::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    std::vector<FileName> command;
    Bool silently;
    Bool continueAfterError;
    Bool executeParallel;

    readConfig(config, "command",            command,            Config::MUSTSET,  "", "");
    readConfig(config, "silently",           silently,           Config::DEFAULT,  "0", "without showing the output.");
    readConfig(config, "continueAfterError", continueAfterError, Config::DEFAULT,  "0", "continue with next command after error, otherwise throw exception");
    readConfig(config, "executeParallel",    executeParallel,    Config::DEFAULT,  "0", "execute several commands in parallel");
    if(isCreateSchema(config)) return;

    // lambda function
    // ---------------
    auto run = [&](UInt i)
    {
      logStatus<<"Run command: \""<<command.at(i)<<"\""<<Log::endl;
      std::vector<std::string> outputs;
      if(!System::exec(command.at(i), outputs))
      {
        if(continueAfterError)
          logWarning<<"Command \""<<command.at(i)<<"\" exited with error"<<Log::endl;
        else
          throw(Exception("Command \""+command.at(i).str()+"\" exited with error"));
      }
      if(!silently)
        for(const auto &output : outputs)
          logInfo<<output<<Log::endl;
    };
    // ---------------

    if(executeParallel)
    {
      Log::GroupPtr groupPtr = Log::group(TRUE, FALSE); // group is freed in the destructor
      Parallel::forEach(command.size(), run, comm, FALSE);
    }
    else if(Parallel::isMaster(comm))
      for(UInt i=0; i<command.size(); i++)
        run(i);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
