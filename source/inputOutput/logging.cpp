/***********************************************/
/**
* @file logging.cpp
*
* @brief Screen and log file output.
*
* @author Torsten Mayer-Guerr
* @author Andreas Kvas
* @date 2007-10-03
*
*/
/***********************************************/

#include "base/import.h"
#include "parallel/parallel.h"
#include "inputOutput/system.h"
#include "logging.h"

/***********************************************/

Log logging;

/***********************************************/

Log::Log()
{
  isLogfile = FALSE;
  printFile = printCout = printCerr = FALSE;
  setSilent(FALSE);
  startTimer();
}

/***********************************************/

Log::~Log()
{
  if(!ss.str().empty())
  {
    std::cerr<<"WARNING: last log line do not end with Log::endl"<<std::endl;
    std::cerr<<"line = '"<<ssFile.str()<<ss.str()<<"'"<<std::endl;
    endLine(ss);
  }
}

/***********************************************/

static Bool isMaster()
{
  return Parallel::isMaster(Parallel::globalCommunicator());
}

/***********************************************/

void Log::setLogFile(const std::string &name)
{
  if(!isMaster()) return;
  if(name.empty())
    return;
  file.open(name, std::ios::out|std::ios::app);
  isLogfile = TRUE;
}

/***********************************************/

void Log::setSilent(Bool silent_)
{
  this->silent = silent_;
}

/***********************************************/

std::ostream &Log::startLine(Type type)
{
  if(!ss.str().empty())
  {
    std::cerr<<"WARNING: last log line do not end with Log::endl"<<std::endl;
    std::cerr<<"line = '"<<ssFile.str()<<ss.str()<<"'"<<std::endl;
    endLine(ss);
  }

  printCerr = ((type==ERROR) || (type==WARNING));
  printCout = (isMaster() && !printCerr);
  printFile = isLogfile && (printCerr || printCout);
  printCout = printCout && !silent;

  if(printFile)
  {
    ssFile<<System::now()%"%y-%m-%d %H:%M:%S"s;
    switch(type)
    {
      case STATUS: ssFile<<" Status  "; break;
      case INFO:   ssFile<<" Info    "; break;
      case WARNING:ssFile<<" WARNING "; break;
      case ERROR:  ssFile<<" ERROR   "; break;
      //default:     ssFile<<"         "; break;
    }
  }

  return ss;
}

/***********************************************/

std::ostream &Log::endLine(std::ostream &stream)
{
  stream<<std::endl;
  if(&stream != &ss)
  {
    std::cerr<<"WARNING: Log::endl used with other ostream than log"<<std::endl;
    return stream;
  }

  if(printCerr)
    std::cerr<<"\033[1;31m"<<ss.str()<<"\033[0m"<<std::flush; // ANSI escape sequence: red and bold
  if(printCout)
    std::cout<<ss.str()<<std::flush;
  if(printFile)
  {
    ssFile<<ss.str();
    file<<ssFile.str()<<std::flush;
  }

  std::cerr<<std::flush;
  std::cout<<std::flush;

  ss.str("");
  ssFile.str("");
  return stream;
}

/***********************************************/
/***********************************************/

void Log::startTimer()
{
  startTime.push(System::now());
}

/***********************************************/

Double Log::seconds() const
{
  return (System::now()-startTime.top()).seconds();
}

/***********************************************/

std::string Log::timeString(Double t)
{
  const Int hour = static_cast<Int>(t)/3600;
  const Int min  = static_cast<Int>(t)/60-60*hour;
  const Int sec  = static_cast<Int>(t)-3600*hour-60*min;
  return hour%"%02i:"s+min%"%02i:"s+sec%"%02i"s;
}

/***********************************************/

void Log::loopTimer(UInt idx, UInt count)
{
  if(isMaster())
  {
    UInt digits = static_cast<UInt>(std::log10(count)) + 1;
    if(count == 0 || (digits >= 4 && ((idx % static_cast<UInt>(std::pow(10, digits-3))) != 0)))
      return;

    const std::string countStr = count%"%i"s;
    Double diff = seconds();
    Double left = (idx!=0) ? diff*(count-idx)/idx : 0;
    if(!silent)
      std::cout<<"\r  "<<std::setw(countStr.size())<<idx+1<<" of "<<countStr<<" (time: "<<timeString(diff)<<", remaining: "<<timeString(left)<<") "<<std::flush;
  }
}

/***********************************************/

void Log::loopTimerEnd(UInt count)
{
  if(isMaster())
  {
    Bool s = silent;
    setSilent(TRUE);
    logStatus<<count<<" loops in time "<<timeString(seconds())<<Log::endl;
    setSilent(s);
    if(!silent)
      std::cout<<"\r  "<<count<<" of "<<count<<" (time: "<<timeString(seconds())<<", remaining: "<<timeString(0.0)<<") "<<std::endl<<std::flush;
  }
  startTime.pop();
}

/***********************************************/
