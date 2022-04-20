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
#include "base/time.h"

/** @addtogroup inputOutputGroup */
/// @{

/***** DEFINES *********************************/

#define logStatus      Log::status()
#define logInfo        Log::info()
#define logWarning     Log::warning()
#define logWarningOnce Log::warningOnce()
#define logError       Log::error()

/***********************************************/

namespace Log
{
  class Group;
  typedef std::shared_ptr<Group> GroupPtr;

  std::function<void(UInt rank, UInt type, const std::string &str)> getReceive();
  void init(UInt rank, UInt size, const std::function<void(UInt type, const std::string &str)> &send);

  GroupPtr group(Bool isMain, Bool silently);  // sub groups of parallel processes
  void setLogFile(const std::string &name);    // set log file for current group (must be called by every process in group)
  void logFilesOnly(Bool enable);              // write only to log file(s) (must be called by every process in group)
  void currentLogFileOnly(Bool enable);        // write only to the current log file in this group (must be called by every process in group)

  std::ostream &status();
  std::ostream &info();
  std::ostream &warningOnce();
  std::ostream &warning();
  std::ostream &error();
  std::ostream &endl(std::ostream &stream);

  class Timer
  {
    Time start;
    UInt count, processCount;
    Bool use;

  public:
    Timer();
    Timer(UInt count, UInt processCount=1, Bool use=TRUE);
    void loopStep(UInt idx);
    void loopEnd() const;
  };
}

/***********************************************/

/// @}

#endif
