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
    LoopPtr loopPtr;
    Bool    continueAfterError;
    Bool    parallelLoops;

    renameDeprecatedConfig(config, "programme", "program", date2time(2020, 6, 3));

    readConfig(config, "loop",               loopPtr,            Config::MUSTSET,  "",  "subprograms are called for every loop");
    readConfig(config, "continueAfterError", continueAfterError, Config::DEFAULT,  "0", "continue with next loop after error, otherwise throw exception");
    readConfig(config, "parallelLoops",      parallelLoops,      Config::DEFAULT,  "0", "parallelize loops instead of programs");
    if(isCreateSchema(config))
    {
      config.xselement("program", "programType", Config::DEFAULT,  Config::UNBOUNDED, "", "");
      return;
    }

    logStatus<<"Run programs"<<Log::endl;
    auto varListOriginal = config.getVarList();
    loopPtr->init(config.getVarList());

    // lambda
    // ------
    auto loopRun = [&]()
    {
      auto varListOld = config.getVarList();
      try
      {
        programRun(config);
      }
      catch(std::exception &e)
      {
        if(!continueAfterError)
        {
          programRemove(config);
          config.getVarList() = varListOriginal;
          throw;
        }
        logError<<e.what()<<"  continue..."<<Log::endl;
      }
      config.getVarList() = varListOld;
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

      logTimerStart;
      for(UInt i=0; !loopPtr->finished(); i++)
      {
        logStatus<<"=== "<<i+1<<". loop ==="<<Log::endl;
        logTimerLoop(i, loopPtr->count());
        loopPtr->setValues(config.getVarList());
        loopRun();
        loopPtr->next(config.getVarList());
      }
      logTimerLoopEnd(loopPtr->count());
    }
    else if(Parallel::isMaster())
    {
      // parallel version: master node
      // -----------------------------
      logTimerStart;
      for(UInt i=0; !loopPtr->finished(); i++)
      {
        logTimerLoop(i, loopPtr->count());
        loopPtr->setValues(config.getVarList());

        UInt process;
        Parallel::receive(process, NULLINDEX); // which process needs work?
        Parallel::send(i, process);            // send new loop number to be computed at process
        loopPtr->next(config.getVarList());
      }
      // send to all processes the end signal (NULLINDEX)
      for(UInt i=1; i<Parallel::size(); i++)
      {
        UInt process;
        Parallel::receive(process, NULLINDEX); // which process needs work?
        Parallel::send(NULLINDEX, process);    // end signal
      }
      Parallel::barrier();
      logTimerLoopEnd(loopPtr->count());
    }
    else
    {
      // clients
      // -------
      Parallel::send(Parallel::myRank(), 0);
      UInt k=0;
      for(;;)
      {
        UInt i;
        Parallel::receive(i, 0);
        if(i==NULLINDEX)
          break;
        for(; k<i; k++) // step to current loop number
          loopPtr->next(config.getVarList());
        loopPtr->setValues(config.getVarList());
        auto comm = Parallel::setDefaultCommunicator(Parallel::selfCommunicator());
        loopRun();
        Parallel::setDefaultCommunicator(comm);
        Parallel::send(Parallel::myRank(), 0);
      }
      Parallel::barrier();
    }

    programRemove(config);
    config.getVarList() = varListOriginal;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
