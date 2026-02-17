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
<<<<<<< HEAD
Determines if a \config{pattern} matches the entire \config{string}.
Supports wildcards * for any number of characters and ? for exactly one character.
If \config{isRegularExpression} is set, \config{pattern} is interpreted as a
regular expression instead. In any case, the \reference{text parser}{general.parser:text}
is applied beforehand.
=======
Determines if a pattern or a regular expression matches the entire string.
Note that, wildcards * and ? are not supported for \config{pattern} when 
\config{isRegularExpression} is 0.
>>>>>>> ec44adb (Add a few words to the doc of the `condition` class)
)";
#endif

/***********************************************/

#include "base/import.h"
#include "base/string.h"
#include "classes/condition/condition.h"
#include <regex>

/***** CLASS ***********************************/

/** @brief Determines if a pattern or a regular expression matches the entire string.
* @ingroup conditionGroup
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
<<<<<<< HEAD
    readConfig(config, "pattern",             pattern,             Config::OPTIONAL, "",  "supports wildcards: * and ?");
=======
    readConfig(config, "pattern",             pattern,             Config::OPTIONAL, "",  "literal '{' in regex should be escaped twice with '##{'");
>>>>>>> ec44adb (Add a few words to the doc of the `condition` class)
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
      p = String::lowerCase(p);
    }

<<<<<<< HEAD
    return std::regex_match(t, ((isRegularExpression) ? std::regex(p) : String::wildcard2regex(p)));
=======
    // Note: because pattern is parsed twice, literal '{' in regex quantifiers must 
    // survive both passes (or avoid braces, e.g. use '^L..$' instead of '^L.{2}$').
    if(isRegularExpression)
      return std::regex_match(t, std::regex(p));
    else
      return (t == p);
>>>>>>> ec44adb (Add a few words to the doc of the `condition` class)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
