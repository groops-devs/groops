/***********************************************/
/**
* @file condition.cpp
*
* @brief Condition.
*
* @author Torsten Mayer-Guerr
* @date 2018-05-18
*
*/
/***********************************************/

#define DOCSTRING_Condition

#include "base/import.h"
#include "config/configRegister.h"
#include "inputOutput/logging.h"
#include "conditionFileExist.h"
#include "conditionCommand.h"
#include "conditionStringContainsPattern.h"
#include "conditionStringMatchPattern.h"
#include "conditionExpression.h"
#include "conditionAnd.h"
#include "conditionOr.h"
#include "conditionNot.h"
#include "condition.h"

/***********************************************/

GROOPS_REGISTER_CLASS(Condition, "conditionType",
                      ConditionFileExist,
                      ConditionCommand,
                      ConditionExpression,
                      ConditionStringContainsPattern,
                      ConditionStringMatchPattern,
                      ConditionAnd,
                      ConditionOr,
                      ConditionNot)

GROOPS_READCONFIG_CLASS(Condition, "conditionType")

/***********************************************/

ConditionPtr Condition::create(Config &config, const std::string &name)
{
  try
  {
    ConditionPtr condition;
    std::string type;
    readConfigChoice(config, name, type, Config::MUSTSET, "", "");

    if(readConfigChoiceElement(config, "fileExist",  type, ""))
      condition = ConditionPtr(new ConditionFileExist(config));
    if(readConfigChoiceElement(config, "command",    type, ""))
      condition = ConditionPtr(new ConditionCommand(config));
    if(readConfigChoiceElement(config, "expression", type, ""))
      condition = ConditionPtr(new ConditionExpression(config));
    if(readConfigChoiceElement(config, "stringContainsPattern", type, ""))
      condition = ConditionPtr(new ConditionStringContainsPattern(config));
    if(readConfigChoiceElement(config, "stringMatchPattern", type, ""))
      condition = ConditionPtr(new ConditionStringMatchPattern(config));
    if(readConfigChoiceElement(config, "and",        type, ""))
      condition = ConditionPtr(new ConditionAnd(config));
    if(readConfigChoiceElement(config, "or",         type, ""))
      condition = ConditionPtr(new ConditionOr(config));
    if(readConfigChoiceElement(config, "not",         type, ""))
      condition = ConditionPtr(new ConditionNot(config));
    endChoice(config);

    return condition;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
