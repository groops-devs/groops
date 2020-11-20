/***********************************************/
/**
* @file loopPrograms.cpp
*
* @brief Runs programs mutiple times.
*
* @author Sebastian Strasser
* @author Torsten Mayer-Guerr
* @date 2017-02-05
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program runs a list of programs in a \configClass{loop}{loopType}.

If \config{continueAfterError}=\verb|yes| and an error occurs, the remaining programs in the current loop
are skipped and iteration continues in the next loop. Otherwise an exception is thrown.
This option is disabled in parallel mode as exceptions cannot be properly caught.

If \config{parallelLoops} is set the loops are distributed to the process nodes, computed in parallel,
and each program is executed with only one node. Otherwise the loops and programs are executed
sequentially but programs are using all nodes.
)";

/***********************************************/

#include "programs/program.h"
#include "classes/loop/loop.h"

/***** CLASS ***********************************/

/** @brief  Runs programs mutiple times.
* @ingroup programsGroup */
class LoopPrograms
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(LoopPrograms, PARALLEL, "Runs programs mutiple times.", System)
GROOPS_RENAMED_PROGRAM(LoopProgramme, LoopPrograms, date2time(2020, 6, 3))

/***********************************************/

void LoopPrograms::run(Config &config)
{
  try
  {
    LoopPtr       loopPtr;
    Bool          continueAfterError;
    Bool          parallelLoops;
    ProgramConfig programs;

    renameDeprecatedConfig(config, "programme", "program", date2time(2020, 6, 3));

    readConfig(config, "loop",               loopPtr,            Config::MUSTSET,  "",  "subprograms are called for every loop");
    readConfig(config, "continueAfterError", continueAfterError, Config::DEFAULT,  "0", "continue with next loop after error, otherwise throw exception");
    readConfig(config, "parallelLoops",      parallelLoops,      Config::DEFAULT,  "0", "parallelize loops instead of programs");
    readConfig(config, "program",            programs,           Config::OPTIONAL, "", "");
    if(isCreateSchema(config)) return;

    logStatus<<"Run programs"<<Log::endl;
    auto varList = config.getVarList();

    // lambda
    // ------
    auto loopRun = [&]()
    {
      try
      {
        auto varListTmp = varList;
        programs.run(varListTmp);
      }
      catch(std::exception &e)
      {
        if(!continueAfterError)
          throw;
        logError<<e.what()<<"  continue..."<<Log::endl;
      }
    };
    // -----------------

    if(!parallelLoops || (Parallel::size() < 3))
    {
      // single process version
      // ----------------------
      if((Parallel::size() > 1) && (continueAfterError))
      {
        if(Parallel::isMaster())
          logWarning<<"continueAfterError does not work on parallel nodes => disabled"<<Log::endl;
        continueAfterError = FALSE;
      }

      UInt iter = 0;
      logTimerStart;
      while(loopPtr->iteration(varList))
      {
        logStatus<<"=== "<<iter+1<<". loop ==="<<Log::endl;
        logTimerLoop(iter++, loopPtr->count());
        loopRun();
      }
      logTimerLoopEnd(loopPtr->count());
      return;
    }

    auto comm = Parallel::setDefaultCommunicator(Parallel::selfCommunicator());
    if(Parallel::isMaster(comm))
    {
      // parallel version: master node
      // -----------------------------
      UInt iter = 0;
      logTimerStart;
      while(loopPtr->iteration(varList))
      {
        logTimerLoop(iter, loopPtr->count());
        UInt process;
        Parallel::receive(process, NULLINDEX, comm); // which process needs work?
        Parallel::send(iter++, process, comm);       // send new loop number to be computed at process
      }
      // send to all processes the end signal (NULLINDEX)
      for(UInt i=1; i<Parallel::size(comm); i++)
      {
        UInt process;
        Parallel::receive(process, NULLINDEX, comm); // which process needs work?
        Parallel::send(NULLINDEX, process, comm);    // end signal
      }
      Parallel::barrier(comm);
      logTimerLoopEnd(loopPtr->count());
    }
    else
    {
      // clients
      // -------
      Parallel::send(Parallel::myRank(comm), 0, comm);
      UInt k=0;
      for(;;)
      {
        UInt i;
        Parallel::receive(i, 0, comm);
        if(i==NULLINDEX)
          break;
        for(; k<=i; k++) // step to current loop number
          loopPtr->iteration(varList);
        loopRun();
        Parallel::send(Parallel::myRank(comm), 0, comm);
      }
      Parallel::barrier(comm);
    }
    Parallel::setDefaultCommunicator(comm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
