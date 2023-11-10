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
  UInt                     index;

public:
  LoopCommandOutput(Config &config);

  UInt count() const override {return strings.size();}
  Bool iteration(VariableList &varList) override;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline LoopCommandOutput::LoopCommandOutput(Config &config)
{
  try
  {
    std::vector<FileName> commands;
    Bool                  silently;

    readConfig(config, "command",            commands,   Config::MUSTSET,   "",            "each output line becomes a loop iteration");
    readConfig(config, "silently",           silently,   Config::DEFAULT,   "0",           "without showing the output.");
    readConfig(config, "variableLoopString", nameString, Config::OPTIONAL, "loopCommand", "name of the variable to be replaced");
    readConfig(config, "variableLoopIndex",  nameIndex,  Config::OPTIONAL, "",            "variable with index of current iteration (starts with zero)");
    readConfig(config, "variableLoopCount",  nameCount,  Config::OPTIONAL, "",            "variable with total number of iterations");
    if(isCreateSchema(config)) return;

    for(const auto &command : commands)
    {
      UInt count = strings.size();
      if(!System::exec(command, strings))
        throw(Exception("Command \""+command.str()+"\" exited with error"));
      if(strings.size() == count)
        logWarningOnce<<"Command \""+command.str()+"\" returned no output"<<Log::endl;
    } // for(i)

    index = 0;
    if(!silently)
      for(const auto &str : strings)
        logInfo<<str<<Log::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Bool LoopCommandOutput::iteration(VariableList &varList)
{
  if(index >= count())
    return FALSE;

  if(!nameString.empty()) varList.setVariable(nameString, strings.at(index));
  if(!nameIndex.empty())  varList.setVariable(nameIndex,  index);
  if(!nameCount.empty())  varList.setVariable(nameCount,  count());

  index++;
  return TRUE;
}

/***********************************************/

#endif
