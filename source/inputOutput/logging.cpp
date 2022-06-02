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

#include <cassert>
#include "base/import.h"
#include "base/string.h"
#include "inputOutput/system.h"
#include "inputOutput/file.h"
#include "logging.h"

/***** CLASS ***********************************/

class Logging
{
public:
  enum Type : UInt {STATUS, INFO, WARNINGONCE, WARNING, ERROR,
                    GROUP, REMOVEGROUP, LOGFILE, LOGFILESONLY, CURRENTLOGFILEONLY};

private:
  struct GroupLocal
  {
    Bool isMain, mustSend;
    GroupLocal(Bool isMain, Bool mustSend) : isMain(isMain), mustSend(mustSend) {}
  };

  // variables at process nodes
  Type                  type;
  std::stringstream     ss;
  std::list<GroupLocal> groupsLocal;
  std::function<void(UInt type, const std::string &str)> send;

  struct Group
  {
    Bool isMain, onScreen, onFile;
    Bool logFilesOnly, currentLogFileOnly;
    std::shared_ptr<OutFile> file;
    Group(Bool isMain, Bool onScreen, Bool onFile) : isMain(isMain), onScreen(onScreen), onFile(onFile), logFilesOnly(FALSE), currentLogFileOnly(FALSE) {}
  };

  // variables at main process
  UInt rank, size;
  Bool newLine;
  std::vector<std::list<Group>> groups; // for each rank/process

  void fatalError(const std::string &str) const;

public:
  Logging();
 ~Logging();

  void init(UInt rank_, UInt size_, const std::function<void(UInt type, const std::string &str)> &send_);

  Log::GroupPtr group(Bool isMain, Bool silently);
  void removeGroup();
  void setLogFile(const std::string &name);
  void logFilesOnly(Bool enable);
  void currentLogFileOnly(Bool enable);

  std::ostream &startLine(Type type);
  std::ostream &endLine(std::ostream &stream);
  void receive(UInt rank, UInt type, const std::string &str);

  friend class Log::Timer;
};

static Logging logging;

/***********************************************/

Logging::Logging() : type(STATUS), groupsLocal({GroupLocal(TRUE, TRUE)}),
                     rank(0), size(1), newLine(FALSE), groups({std::list<Group>({Group(TRUE, TRUE, TRUE)})})
{
  send = std::bind(&Logging::receive, this, 0, std::placeholders::_1, std::placeholders::_2);
}

/***********************************************/

Logging::~Logging()
{
  if(!ss.str().empty())
  {
    endLine(ss);
    std::cerr<<"WARNING: last log line does not end with with Log::endl"<<std::endl;
    std::cerr<<"line = '"<<ss.str()<<"'"<<std::endl;
  }
}

/***********************************************/

void Logging::init(UInt rank_, UInt size_, const std::function<void(UInt type, const std::string &str)> &send_)
{
  try
  {
    rank = rank_;
    size = size_;
    send = send_;

    groupsLocal.clear();
    groupsLocal.emplace_front((rank == 0), TRUE);

    groups.clear();
    if(rank == 0)
      groups.resize(size, std::list<Group>({Group((rank == 0), TRUE, TRUE)}));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

class Log::Group
{
public:
  Group() {}
 ~Group() {logging.removeGroup();}
};

/***********************************************/

Log::GroupPtr Logging::group(Bool isMain, Bool silently)
{
  try
  {
    assert(groupsLocal.size());
    send(GROUP, (isMain ? "1"s : "0"s) + (!silently ? "1"s : "0"s));
    const Bool mustSendBefore = groupsLocal.front().mustSend;
    groupsLocal.emplace_front(isMain, !silently && mustSendBefore);
    return std::make_shared<Log::Group>();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Logging::removeGroup()
{
  try
  {
    assert(groupsLocal.size());
    send(REMOVEGROUP, "");
    groupsLocal.pop_front();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Logging::setLogFile(const std::string &name)
{
  try
  {
    if(name.empty())
      return;
    assert(groupsLocal.size());
    send(LOGFILE, name);
    groupsLocal.front().mustSend = TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Logging::logFilesOnly(Bool enable)
{
  try
  {
    send(LOGFILESONLY, enable ? "1" : "0");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Logging::currentLogFileOnly(Bool enable)
{
  try
  {
    send(CURRENTLOGFILEONLY, enable ? "1" : "0");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

std::ostream &Logging::startLine(Type type_)
{
  try
  {
    if(!ss.str().empty())
    {
      endLine(ss);
      std::cerr<<"WARNING: last log line does not end with with Log::endl"<<std::endl;
      std::cerr<<"line = '"<<ss.str()<<"'"<<std::endl;
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
    assert(groupsLocal.size());
    if((groupsLocal.front().isMain && (groupsLocal.front().mustSend || (type == WARNINGONCE))) || (type == WARNING) || (type == ERROR))
      for(const std::string &str :  String::split(ss.str(), '\n'))
        send(type, str);

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
void Logging::receive(UInt rank, UInt type, const std::string &str)
{
  try
  {
    assert(groups.size() && groups.at(rank).size());

    if(type == GROUP)
    {
      const Bool onScreen = groups.at(rank).front().onScreen;
      groups.at(rank).emplace_front(str.at(0) == '1', (str.at(1) == '1') && onScreen, (str.at(1) == '1'));
      return;
    }

    if(type == REMOVEGROUP)
    {
      groups.at(rank).pop_front();
      return;
    }

    if(type == LOGFILE)
    {
      // already open?
      for(auto &groupsRank : groups)
        for(auto &group : groupsRank)
          if(group.file && (group.file->fileName().str() == str))
          {
            groups.at(rank).front().file = group.file;
            return;
          }
      groups.at(rank).front().file = std::make_shared<OutFile>(str);
      return;
    }

    if(type == LOGFILESONLY)
    {
      groups.at(rank).front().logFilesOnly = (str == "1");
      return;
    }

    if(type == CURRENTLOGFILEONLY)
    {
      groups.at(rank).front().currentLogFileOnly = (str == "1");
      return;
    }

    // screen
    // ------
    if(!groups.at(rank).front().logFilesOnly && !groups.at(rank).front().currentLogFileOnly)
      if((groups.at(rank).front().isMain && (groups.at(rank).front().onScreen || (type == WARNINGONCE))) || (type == WARNING) || (type == ERROR))
      {
        if(newLine)
          std::cout<<std::endl<<std::flush;
        newLine = FALSE;
        if((type == STATUS) || (type == INFO))
          std::cout<<(rank ? rank%"(process%3i) "s : ""s)<<str<<std::endl<<std::flush;
        else
#ifdef _WIN32
          std::cerr<<(rank ? rank%"(process%3i) "s : ""s)<<str<<std::endl<<std::flush;
#else
          std::cerr<<"\033[1;31m"<<(rank ? rank%"(process%3i) "s : ""s)<<str<<"\033[0m"<<std::endl<<std::flush; // ANSI escape sequence: red and bold
#endif
      }

    // log files
    // ---------
    for(auto &group : groups.at(rank))
    {
      if(group.file)
        switch(type)
        {
          case STATUS:      *(group.file)<<System::now()%"%y-%m-%d %H:%M:%S"s<<" Status  "<<(group.isMain ? ""s : rank%"(process%3i) "s)<<str<<std::endl; break;
          case INFO:        *(group.file)<<System::now()%"%y-%m-%d %H:%M:%S"s<<" Info    "<<(group.isMain ? ""s : rank%"(process%3i) "s)<<str<<std::endl; break;
          case WARNINGONCE: *(group.file)<<System::now()%"%y-%m-%d %H:%M:%S"s<<" WARNING "<<(group.isMain ? ""s : rank%"(process%3i) "s)<<str<<std::endl; break;
          case WARNING:     *(group.file)<<System::now()%"%y-%m-%d %H:%M:%S"s<<" WARNING "<<(group.isMain ? ""s : rank%"(process%3i) "s)<<str<<std::endl; break;
          case ERROR:       *(group.file)<<System::now()%"%y-%m-%d %H:%M:%S"s<<" ERROR   "<<(group.isMain ? ""s : rank%"(process%3i) "s)<<str<<std::endl; break;
        }
      // send to logfile in higher level?
      if(group.currentLogFileOnly)
        break;
      if(!(group.isMain && (group.onFile || (type == WARNINGONCE))) && (type != WARNING) && (type != ERROR))
        break;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

Log::Timer::Timer() : start(System::now()), count(0), processCount(1), use(TRUE) {}

/***********************************************/

Log::Timer::Timer(UInt count, UInt processCount, Bool use) : start(System::now()), count(count), processCount(processCount), use(use) {}

/***********************************************/

void Log::Timer::loopStep(UInt idx)
{
  try
  {
    assert((logging.rank > 0) || (logging.groups.size() && logging.groups.at(0).size()));
    if(!use || (count == 0) || (logging.rank != 0) || !logging.groups.at(0).front().onScreen)
      return;

    if(idx >= 5000)
    {
      const UInt digits = static_cast<UInt>(std::log10(idx+1)) + 1;
      if((idx+1) % static_cast<UInt>(std::pow(10, digits-3)))
        return;
    }

    count = std::max(count, idx+1);
    const Double diff    = (System::now()-start).mjd();
    const Double perStep = (idx >= processCount) ? diff/(idx/processCount) : 0.;
    const Double left    = ((count+processCount-1)/processCount - idx/processCount) * perStep;

    const std::string countStr = count%"%i"s;
    std::cout<<"\r  "<<std::setw(countStr.size())<<idx+1<<" of "<<countStr<<" (time: "<<diff%"%H:%M:%S, remaining: "s<<left%"%H:%M:%S) "s<<std::flush;
    logging.newLine = TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Log::Timer::loopEnd() const
{
  try
  {
    if(!use)
      return;

    assert((logging.rank > 0) || (logging.groups.size() && logging.groups.at(0).size()));
    if((logging.rank == 0) && logging.groups.at(0).front().onScreen)
    {
      std::cout<<"\r  "<<count<<" of "<<count<<" (time: "<<(System::now()-start)%"%H:%M:%S, remaining: "s<<0.%"%H:%M:%S) "s<<std::endl<<std::flush;
      logging.newLine = FALSE;
    }
    // to log file(s)
    logFilesOnly(TRUE);
    logStatus<<count<<" loops in time "<<(System::now()-start)%"%H:%M:%S"s<<Log::endl;
    logFilesOnly(FALSE);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

std::function<void(UInt rank, UInt type, const std::string &str)> Log::getReceive() {return std::bind(&Logging::receive, &logging, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);}
void Log::init(UInt rank, UInt size, const std::function<void(UInt type, const std::string &str)> &send) {logging.init(rank, size, send);}
Log::GroupPtr Log::group(Bool isMain, Bool silently) {return logging.group(isMain, silently);}
void Log::setLogFile(const std::string &name)        {logging.setLogFile(name);}
void Log::logFilesOnly(Bool enable)                  {logging.logFilesOnly(enable);}
void Log::currentLogFileOnly(Bool enable)            {logging.currentLogFileOnly(enable);}
std::ostream &Log::status()                          {return logging.startLine(Logging::STATUS);}
std::ostream &Log::info()                            {return logging.startLine(Logging::INFO);}
std::ostream &Log::warningOnce()                     {return logging.startLine(Logging::WARNINGONCE);}
std::ostream &Log::warning()                         {return logging.startLine(Logging::WARNING);}
std::ostream &Log::error()                           {return logging.startLine(Logging::ERROR);}
std::ostream &Log::endl(std::ostream &stream)        {return logging.endLine(stream);}

/***********************************************/
