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
#include "base/string.h"
#include "inputOutput/system.h"
#include "inputOutput/file.h"
#include "logging.h"

/***** CLASS ***********************************/

class Logging
{
public:
  enum Type : UInt {STATUS, INFO, WARNINGONCE, WARNING, ERROR};

  Type              type;
  UInt              rank;
  Bool              enabled, silent, newLine, isLogfile;
  std::stringstream ss;
  std::function<void(UInt type, const std::string &str)> send;
  OutFile           file;

  Logging();
 ~Logging();

  void setRank(UInt rank_)       {rank = rank_;}
  Bool enableOutput(Bool enable) {std::swap(enable, enabled);  return enable;}
  void setSilent(Bool silent_)   {silent = silent_;}
  void setLogFile(const std::string &name);

  std::ostream &startLine(Type type);
  std::ostream &endLine(std::ostream &stream);
  void receive(UInt type, const std::string &str);

  // Timer
  std::stack<Time> startTime;
  void   startTimer();
  void   loopTimer(UInt idx, UInt count, UInt processCount);
  void   loopTimerEnd(UInt count);
};

static Logging logging;

/***********************************************/

Logging::Logging() : type(STATUS), rank(0), enabled(TRUE), silent(FALSE), newLine(FALSE)
{
  send = std::bind(&Logging::receive, this, std::placeholders::_1, std::placeholders::_2);
  startTimer();
}

/***********************************************/

Logging::~Logging()
{
  if(!ss.str().empty())
  {
    std::cerr<<"WARNING: last log line does not end with with Log::endl"<<std::endl;
    std::cerr<<"line = '"<<ss.str()<<"'"<<std::endl;
    endLine(ss);
  }
}

/***********************************************/

void Logging::setLogFile(const std::string &name)
{
  if(rank != 0) return;
  if(name.empty())
    return;
  file.open(name);
  isLogfile = TRUE;
}

/***********************************************/

std::ostream &Logging::startLine(Type type_)
{
  try
  {
    if(!ss.str().empty())
    {
      std::cerr<<"WARNING: last log line does not end with with Log::endl"<<std::endl;
      std::cerr<<"line = '"<<ss.str()<<"'"<<std::endl;
      endLine(ss);
    }

    type = type_;
    return ss;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::ostream &Logging::endLine(std::ostream &stream)
{
  try
  {
    if(&stream != &ss)
      throw(Exception("Log::endl used with other ostream than log"));

    // send log line to main process
    if(enabled || (type == WARNING) || (type == ERROR))
      for(const std::string &str :  String::split(ss.str(), '\n'))
      {
        if(rank == 0)
          receive(type, str);
        else
          send(type, rank%"%4i. process: "s+str);
      }

    ss.str("");
    return stream;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

// received log line at main node
void Logging::receive(UInt type, const std::string &str)
{
  try
  {
    if(!silent || (type == WARNINGONCE) || (type == WARNING) || (type == ERROR))
    {
      if(newLine)
        std::cout<<std::endl<<std::flush;
      newLine = FALSE;
      if((type == STATUS) || (type == INFO))
        std::cout<<str<<std::endl<<std::flush;
      else
        std::cerr<<"\033[1;31m"<<str<<"\033[0m"<<std::endl<<std::flush; // ANSI escape sequence: red and bold
    }

    if(isLogfile)
    {
      file<<System::now()%"%y-%m-%d %H:%M:%S"s;
      switch(type)
      {
        case STATUS:      file<<" Status  "; break;
        case INFO:        file<<" Info    "; break;
        case WARNINGONCE: file<<" WARNING "; break;
        case WARNING:     file<<" WARNING "; break;
        case ERROR:       file<<" ERROR   "; break;
      }
      file<<str<<std::endl;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void Logging::startTimer()
{
  try
  {
    startTime.push(System::now());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Logging::loopTimer(UInt idx, UInt count, UInt processCount)
{
  try
  {
    if((count == 0) || (rank != 0) || !enabled || silent)
      return;

    if(idx >= 5000)
    {
      const UInt digits = static_cast<UInt>(std::log10(idx+1)) + 1;
      if((idx+1) % static_cast<UInt>(std::pow(10, digits-3)))
        return;
    }

    const Double diff    = (System::now()-startTime.top()).mjd();
    const Double perStep = (idx >= processCount) ? diff/(idx/processCount) : 0.;
    const Double left    = (count-idx+processCount-1)/processCount * perStep;

    const std::string countStr = count%"%i"s;
    std::cout<<"\r  "<<std::setw(countStr.size())<<idx+1<<" of "<<countStr<<" (time: "<<diff%"%H:%M:%S, remaining: "s<<left%"%H:%M:%S) "s<<std::flush;
    newLine = TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Logging::loopTimerEnd(UInt count)
{
  try
  {
    Bool s = silent;
    if((rank == 0) && enabled && !silent)
    {
      std::cout<<"\r  "<<count<<" of "<<count<<" (time: "<<(System::now()-startTime.top())%"%H:%M:%S, remaining: "s<<0.%"%H:%M:%S) "s<<std::endl<<std::flush;
      newLine = FALSE;
      silent  = TRUE;
    }
    logStatus<<count<<" loops in time "<<(System::now()-startTime.top())%"%H:%M:%S"s<<Log::endl;
    setSilent(s);
    startTime.pop();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

std::function<void(UInt type, const std::string &str)> Log::getReceive() {return std::bind(&Logging::receive, &logging, std::placeholders::_1, std::placeholders::_2);}
void Log::setSend(std::function<void(UInt type, const std::string &str)> send) {logging.send = send;}
void Log::setRank(UInt rank)                                 {logging.setRank(rank);}
Bool Log::enableOutput(Bool enable)                          {return logging.enableOutput(enable);}
void Log::setSilent(Bool silent)                             {logging.setSilent(silent);}
void Log::setLogFile(const std::string &name)                {logging.setLogFile(name);}
void Log::startTimer()                                       {logging.startTimer();}
void Log::loopTimer(UInt idx, UInt count, UInt processCount) {logging.loopTimer(idx, count, processCount);}
void Log::loopTimerEnd(UInt count)                           {logging.loopTimerEnd(count);}
std::ostream &Log::status()                                  {return logging.startLine(Logging::STATUS);}
std::ostream &Log::info()                                    {return logging.startLine(Logging::INFO);}
std::ostream &Log::warningOnce()                             {return logging.startLine(Logging::WARNINGONCE);}
std::ostream &Log::warning()                                 {return logging.startLine(Logging::WARNING);}
std::ostream &Log::error()                                   {return logging.startLine(Logging::ERROR);}
std::ostream &Log::endl(std::ostream &stream)                {return logging.endLine(stream);}

/***********************************************/
