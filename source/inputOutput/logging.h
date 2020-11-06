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
* Next to logStatus there are also logInfo, logWarning, logError and logDebug(level).
* logDebug messages are only output if setDebugLevel(level) is at least level.
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
#include "inputOutput/file.h"

/** @addtogroup inputOutputGroup */
/// @{

/***** CLASS ***********************************/

class  Log;
extern Log logging;

/**
* @brief Screen and log file output.
*
* Calling logStatus<<"Text"<<Log::endl;
* outputs "Text" on the screen
* and "2007-10-03 11:49:35 Status Text" in the log file.
* Only outputs from the main process are printed.
* Next to logStatus there are also logInfo, logWarning, logError and logDebug(level).
* logDebug messages are only output if setDebugLevel(level) is at least level.
* IMPORTANT: Each output must end with Log::endl.
*/
class Log
{
public:
  enum Type {STATUS, INFO, WARNING, ERROR};

  Log();
 ~Log();

  Log(const Log &) = delete;
  Log &operator=(const Log &) = delete;

  void setLogFile(const std::string &name);
  void setSilent(Bool silent=TRUE);

  // Timer
  void   startTimer();
  void   loopTimer(UInt idx, UInt count);
  void   loopTimerEnd(UInt count);
  Double seconds() const;
  static std::string timeString(Double t);

  // Internal functions
  std::ostream &startLine(Type type=STATUS);
  std::ostream &endLine(std::ostream &stream);
  static std::ostream &endl(std::ostream &stream) {return logging.endLine(stream);}

private:
  std::stringstream ss, ssFile;
  OutFile           file;
  Bool              isLogfile;
  Bool              silent;
  Bool              printFile, printCout, printCerr;
  std::stack<Time>  startTime; // Timer
};

/***********************************************/

#define logStatus       logging.startLine(Log::STATUS)
#define logInfo         logging.startLine(Log::INFO)
#define logWarning      logging.startLine(Log::WARNING)
#define logError        logging.startLine(Log::ERROR)

#define logTimerStart           logging.startTimer();
#define logTimerLoop(idx,count) logging.loopTimer(idx, count);
#define logTimerLoopEnd(count)  logging.loopTimerEnd(count);

/***********************************************/

/// @}

#endif
