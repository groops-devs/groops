/***********************************************/
/**
* @file logging.h
*
* @brief Screen and log file output.
*
* Calling logStatus<<"Text"<<Log::endl;
* outputs "Text" on the screen
* and "2007-10-03 11:49:35 Status Text" in the log file.
* Only outputs from the main process are printed.
* Next to logStatus there are also logInfo, logWarning, logError.
* IMPORTANT: Each output must end with Log::endl.
*
* @author Torsten Mayer-Guerr
* @author Andreas Kvas
* @date 2007-10-03
*
*/
/***********************************************/

#ifndef __GROOPS_LOGGING__
#define __GROOPS_LOGGING__

#include "base/importStd.h"

/** @addtogroup inputOutputGroup */
/// @{

/***** DEFINES *********************************/

#define logStatus               Log::status()
#define logInfo                 Log::info()
#define logWarning              Log::warning()
#define logWarningOnce          Log::warningOnce()
#define logError                Log::error()
#define logTimerStart           Log::startTimer();
#define logTimerLoop(idx,count) Log::loopTimer(idx, count);
#define logTimerLoopEnd(count)  Log::loopTimerEnd(count);

/***********************************************/

namespace Log
{
  std::function<void(UInt type, const std::string &str)> getRecieve();
  void setSend(std::function<void(UInt type, const std::string &str)> send);
  void setRank(UInt rank);
  Bool enableOutput(Bool enable);
  void setSilent(Bool silent);
  void setLogFile(const std::string &name);

  std::ostream &status();
  std::ostream &info();
  std::ostream &warningOnce();
  std::ostream &warning();
  std::ostream &error();
  std::ostream &endl(std::ostream &stream);

  void startTimer();
  void loopTimer(UInt idx, UInt count, UInt processCount=1);
  void loopTimerEnd(UInt count);
}

/***********************************************/

/// @}

#endif
