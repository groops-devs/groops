/***********************************************/
/**
* @file loopPrograms.cpp
*
* @brief Runs programs multiple times.
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

If \config{continueAfterError}=\verb|yes| and an error occurs, the remaining programs in the current iteration
are skipped and the loop continues with the next iteration. Otherwise an exception is thrown.

If this program is executed on multiple processing nodes, the iterations can be computed in parallel,
see \reference{parallelization}{general.parallelization}. The first process serves as load balancer
and the other processes are assigned to iterations according to \config{processCountPerIteration}.
For example, running a loop containing three iterations on 13 processes with \config{processCountPerIteration}=\verb|4|,
runs the three iterations in parallel, with each iteration being assigned four processes.
With \config{parallelLog}=\verb|yes| all processes write output to screen and the log file.
As the output can be quite confusing in this case, running \program{GroupPrograms} with an extra \config{outputfileLog}
for each iteration (use the loop variables for the name of the log files) might be helpful.
)";

/***********************************************/

#include "programs/program.h"
#include "classes/loop/loop.h"

/***** CLASS ***********************************/

/** @brief  Runs programs multiple times.
* @ingroup programsGroup */
class LoopPrograms
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(LoopPrograms, PARALLEL, "Runs programs multiple times.", System)
GROOPS_RENAMED_PROGRAM(LoopProgramme, LoopPrograms, date2time(2020, 6, 3))

/***********************************************/

void LoopPrograms::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    LoopPtr       loopPtr;
    Bool          continueAfterError;
    UInt          processCount;
    Bool          parallelLog;
    ProgramConfig programs;

    renameDeprecatedConfig(config, "programme",     "program", date2time(2020, 6, 3));
    renameDeprecatedConfig(config, "parallelLoops", "processCountPerLoopStep", date2time(2020, 12, 28));
    renameDeprecatedConfig(config, "processCountPerLoopStep", "processCountPerIteration", date2time(2020, 1, 14));

    readConfig(config, "loop",                     loopPtr,            Config::MUSTSET,  "",  "subprograms are called for every iteration");
    readConfig(config, "continueAfterError",       continueAfterError, Config::DEFAULT,  "0", "continue with next iteration after error, otherwise throw exception");
    readConfig(config, "processCountPerIteration", processCount,       Config::DEFAULT,  "0", "0: use all processes for each iteration");
    readConfig(config, "parallelLog",              parallelLog,        Config::DEFAULT,  "1", "write to screen/log file from all processing nodes in parallelized loops");
    readConfig(config, "program",                  programs,           Config::OPTIONAL, "",  "");
    if(isCreateSchema(config)) return;

    VariableList varList;

    // Every process executes every iteration
    // --------------------------------------
    if((processCount == 0) || (processCount+1 >= Parallel::size(comm)) || (Parallel::size(comm) < 3))
    {
      UInt iter = 0;
      if(loopPtr->iteration(varList))
      {
        Log::Timer timer(loopPtr->count());
        do
        {
          logStatus<<"=== "<<iter+1<<". loop ==="<<Log::endl;
          timer.loopStep(iter++);
          try
          {
            Parallel::broadCastExceptions(comm, [&](Parallel::CommunicatorPtr comm)
            {
              auto varListTmp = varList;
              programs.run(varListTmp, comm);
            });
          }
          catch(std::exception &e)
          {
            if(!continueAfterError || Parallel::isExternal(e))
              throw;
            if(Parallel::isMaster(comm))
              logError<<e.what()<<"  continue..."<<Log::endl;
          }
        }
        while(loopPtr->iteration(varList));
        timer.loopEnd();
      }
      return;
    }

    // Iterations in parallel
    // ----------------------
    processCount = std::min(processCount, Parallel::size(comm)-1);
    UInt rank = Parallel::myRank(comm);
    auto commLocal = Parallel::splitCommunicator((rank > 0) ? (rank-1)/processCount : NULLINDEX, rank, comm); // processes of an iteration
    auto commLoop  = Parallel::splitCommunicator(((rank == 0) || ((rank-1)%processCount == 0)) ?  0 : NULLINDEX, rank, comm); // 'main' processes of all iterations
    Log::GroupPtr groupPtr = Log::group(((rank == 0) || ((rank-1)%processCount == 0)), (rank>0) && !parallelLog); // group is freed in the destructor

    if(commLoop && Parallel::isMaster(commLoop))
    {
      // parallel version: main node
      // ---------------------------
      UInt iter = 0;
      if(loopPtr->iteration(varList))
      {
        Log::Timer timer(loopPtr->count(), Parallel::size(commLoop)-1, TRUE);
        do
        {
          timer.loopStep(iter);
          UInt process;
          Parallel::receive(process, NULLINDEX, commLoop); // which process needs work?
          Parallel::send(iter++, process, commLoop);       // send new loop number to be computed at process
        }
        while(loopPtr->iteration(varList));
        // send to all processes the end signal (NULLINDEX)
        for(UInt i=1; i<Parallel::size(commLoop); i++)
        {
          UInt process;
          Parallel::receive(process, NULLINDEX, commLoop); // which process needs work?
          Parallel::send(NULLINDEX, process, commLoop);    // end signal
        }
        Parallel::barrier(comm);
        timer.loopEnd();
      }
    }
    else
    {
      // clients
      // -------
      UInt k=0;
      for(;;)
      {
        UInt i;
        if(Parallel::isMaster(commLocal))
        {
          Parallel::send(Parallel::myRank(commLoop), 0, commLoop);
          Parallel::receive(i, 0, commLoop);
        }
        Parallel::broadCast(i, 0, commLocal);
        if(i == NULLINDEX) // end signal?
          break;
        for(; k<=i; k++) // step to current loop number
          loopPtr->iteration(varList);

        try
        {
          Parallel::broadCastExceptions(commLocal, [&](Parallel::CommunicatorPtr commLocal)
          {
            auto varListTmp = varList;
            programs.run(varListTmp, commLocal);
          });
        }
        catch(std::exception &e)
        {
          if(!continueAfterError || Parallel::isExternal(e))
            throw;
          if(Parallel::isMaster(commLocal))
            logError<<e.what()<<"  continue..."<<Log::endl;
        }
      }
      Parallel::barrier(comm);
    } // clients
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
