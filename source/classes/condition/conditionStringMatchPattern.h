/***********************************************/
/**
* @file conditionStringMatchPattern.h
*
* @brief Determines if a pattern or a regular expression matches the entire string.
*
* @author Torsten Mayer-Guerr
* @date 2019-03-17
*
*/
/***********************************************/

#ifndef __GROOPS_CONDITIONSTRINGMATCHPATTERN__
#define __GROOPS_CONDITIONSTRINGMATCHPATTERN__

// Latex documentation
#ifdef DOCSTRING_Condition
static const char *docstringConditionStringMatchPattern = R"(
\subsection{StringMatchPattern}
Determines if a pattern or a regular expression matches the entire string.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "base/string.h"
#include "classes/condition/condition.h"
#include <regex>

/***** CLASS ***********************************/

/** @brief Determines if a pattern or a regular expression matches the entire string.
* @ingroup ConditionGroup
* @see Condition */
class ConditionStringMatchPattern : public Condition
{
  FileName text, pattern;
  Bool     isRegularExpression, caseSensitive;

public:
  ConditionStringMatchPattern(Config &config);

  Bool condition(const VariableList &varList) const;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline ConditionStringMatchPattern::ConditionStringMatchPattern(Config &config)
{
  try
  {
    readConfig(config, "string",              text,                Config::OPTIONAL, "",  "should contain a {variable}");
    readConfig(config, "pattern",             pattern,             Config::OPTIONAL, "",  "");
    readConfig(config, "isRegularExpression", isRegularExpression, Config::DEFAULT,  "0", "pattern is  a regular expression");
    readConfig(config, "caseSensitive",       caseSensitive,       Config::DEFAULT,  "1", "treat lower and upper case as distinct");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Bool ConditionStringMatchPattern::condition(const VariableList &varList) const
{
  try
  {
    std::string t = text(varList).str();
    std::string p = pattern(varList).str();
    if(!caseSensitive)
    {
      t = String::lowerCase(t);
      p = String::lowerCase(t);
    }

    if(isRegularExpression)
      return std::regex_match(t, std::regex(p));
    else
      return (t == p);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
