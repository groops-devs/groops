/***********************************************/
/**
* @file conditionStringContainsPattern.h
*
* @brief Determines if there is a match between a pattern or a regular expression and some subsequence in a string.
*
* @author Torsten Mayer-Guerr
* @date 2019-03-17
*
*/
/***********************************************/

#ifndef __GROOPS_CONDITIONSTRINGCONTAINSPATTERN__
#define __GROOPS_CONDITIONSTRINGCONTAINSPATTERN__

// Latex documentation
#ifdef DOCSTRING_Condition
static const char *docstringConditionStringContainsPattern = R"(
\subsection{StringContainsPattern}
Determines if there is a match between a pattern or a regular expression and some subsequence in a string.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "base/string.h"
#include "classes/condition/condition.h"
#include <regex>

/***** CLASS ***********************************/

/** @brief Determines if there is a match between a pattern or a regular expression and some subsequence in a string.
* @ingroup conditionGroup
* @see Condition */
class ConditionStringContainsPattern : public Condition
{
  FileName text, pattern;
  Bool     isRegularExpression, caseSensitive;

public:
  ConditionStringContainsPattern(Config &config);

  Bool condition(const VariableList &varList) const;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline ConditionStringContainsPattern::ConditionStringContainsPattern(Config &config)
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

inline Bool ConditionStringContainsPattern::condition(const VariableList &varList) const
{
  try
  {
    std::string t = text(varList).str();
    std::string p = pattern(varList).str();
    if(!caseSensitive)
    {
      t = String::lowerCase(t);
      p = String::lowerCase(p);
    }

    if(isRegularExpression)
      return std::regex_search(t, std::regex(p));
    else
      return (t.find(p) != std::string::npos);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
