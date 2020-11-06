/***********************************************/
/**
* @file loopCommandOutput.h
*
* @brief Loop over lines of command output.
*
* @author Sebastian Strasser
* @date 2017-02-10
*
*/
/***********************************************/

#ifndef __GROOPS_LOOPCOMMANDOUTPUT__
#define __GROOPS_LOOPCOMMANDOUTPUT__

// Latex documentation
#ifdef DOCSTRING_Loop
static const char *docstringLoopCommandOutput = R"(
\subsection{CommandOutput}\label{loopType:commandOutput}
Loop over lines of command output.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "parallel/parallel.h"
#include "inputOutput/logging.h"
#include "inputOutput/system.h"
#include "classes/loop/loop.h"

/***** CLASS ***********************************/

/** @brief Loop over lines of command output.
* @ingroup LoopGroup
* @see Loop */
class LoopCommandOutput : public Loop
{
  std::string              nameString, nameIndex, nameCount;
  std::vector<std::string> strings;
  std::vector<FileName>    command;
  Bool                     silently;

public:
  LoopCommandOutput(Config &config);

  UInt count() const override;
  void init(VariableList &varList) override;
  void setValues(VariableList &varList) override;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline LoopCommandOutput::LoopCommandOutput(Config &config) : Loop()
{
  try
  {
    readConfig(config, "command",            command,    Config::MUSTSET,   "",            "each output line becomes a loop iteration");
    readConfig(config, "silently",           silently,   Config::DEFAULT,   "0",           "without showing the output.");
    readConfig(config, "variableLoopString", nameString, Config::OPTIONAL, "loopCommand", "name of the variable to be replaced");
    readConfig(config, "variableLoopIndex",  nameIndex,  Config::OPTIONAL, "",            "variable with index of current iteration (starts with zero)");
    readConfig(config, "variableLoopCount",  nameCount,  Config::OPTIONAL, "",            "variable with total number of iterations");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline UInt LoopCommandOutput::count() const
{
  return strings.size();
}

/***********************************************/

inline void LoopCommandOutput::init(VariableList &varList)
{
  try
  {
    strings.clear();
    if(Parallel::isMaster())
    {
      UInt count = 0;
      for(UInt i=0; i<command.size(); i++)
      {
        if(!System::exec(command.at(i)(varList), strings))
          throw(Exception("Command \""+command.at(i)(varList).str()+"\" exited with error"));
        if(!silently)
          for(const auto &str : strings)
            logInfo<<str<<Log::endl;
        if(strings.size() == count)
          logWarning<<"Command \""+command.at(i)(varList).str()+"\" returned no output"<<Log::endl;
        count = strings.size();
      } // for(i)
    }
    Parallel::broadCast(strings);

    if(index == NULLINDEX)
    {
      if(!nameString.empty()) addVariable(nameString, varList);
      if(!nameIndex.empty())  addVariable(nameIndex,  varList);
      if(!nameCount.empty())  addVariable(nameCount,  varList);
    }
    index = 0;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void LoopCommandOutput::setValues(VariableList &varList)
{
  if(!nameString.empty()) varList[nameString]->setValue(strings.at(index));
  if(!nameIndex.empty())  varList[nameIndex]->setValue(index);
  if(!nameCount.empty())  varList[nameCount]->setValue(count());
}

/***********************************************/

#endif
