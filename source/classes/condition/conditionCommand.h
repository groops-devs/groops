/***********************************************/
/**
* @file conditionCommand.h
*
* @brief Execute command and check success.
*
* @author Torsten Mayer-Guerr
* @date 2018-08-18
*
*/
/***********************************************/

#ifndef __GROOPS_CONDITIONCOMMAND__
#define __GROOPS_CONDITIONCOMMAND__

// Latex documentation
#ifdef DOCSTRING_Condition
static const char *docstringConditionCommand = R"(
\subsection{Command}
Execute command and check success.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "parallel/parallel.h"
#include "parser/stringParser.h"
#include "inputOutput/system.h"
#include "classes/condition/condition.h"

/***** CLASS ***********************************/

/** @brief Execute command and check success.
* @ingroup ConditionGroup
* @see Condition */
class ConditionCommand : public Condition
{
  FileName command;
  Bool     silently;

public:
  ConditionCommand(Config &config);

  Bool condition(const VariableList &varList) const;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline ConditionCommand::ConditionCommand(Config &config)
{
  try
  {
    readConfig(config, "command",  command,  Config::MUSTSET,  "",  "");
    readConfig(config, "silently", silently, Config::DEFAULT,  "0", "without showing the output.");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Bool ConditionCommand::condition(const VariableList &varList) const
{
  try
  {
    std::vector<std::string> outputs;
    Bool result = System::exec(StringParser::parse(command, varList), outputs);
    if(!silently)
      for(const auto &output : outputs)
        logInfo<<output<<Log::endl;
  return result;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
