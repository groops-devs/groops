/***********************************************/
/**
* @file string.cpp
*
* @brief Miscellaneous string functions.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2019-09-04
*
*/
/***********************************************/

#include "base/importStd.h"
#include "base/string.h"
#include <regex>

/***********************************************/

std::string String::upperCase(const std::string &str)
{
  std::string result = str;
  std::transform(result.begin(), result.end(), result.begin(), [](auto c){return std::toupper(c);});
  return result;
}

/***********************************************/

std::string String::lowerCase(const std::string &str)
{
  std::string result = str;
  std::transform(result.begin(), result.end(), result.begin(), [](auto c){return std::tolower(c);});
  return result;
}

/***********************************************/

std::string String::trim(const std::string &str)
{
  try
  {
    auto start = str.find_first_not_of(" \t\r");
    if(start == std::string::npos)
      return "";
    auto end = str.find_last_not_of(" \t\r");
    return str.substr(start, end-start+1);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

std::string String::trimLeft(const std::string &str)
{
  try
  {
    auto start = str.find_first_not_of(" \t\r");
    if(start == std::string::npos)
      return "";
    return str.substr(start);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

std::string String::trimRight(const std::string &str)
{
  try
  {
    auto end = str.find_last_not_of(" \t\r");
    if(end == std::string::npos)
      return "";
    return str.substr(0, end+1);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

Double String::toDouble(const std::string &str)
{
  try
  {
    if(std::all_of(str.begin(), str.end(), ::isspace))
      return 0.;

    std::string str2 = str;
    auto dpos = str.find_first_of("Dd");
    if(dpos != std::string::npos)
      str2[dpos] = 'e';

    return std::stod(str2);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("cannot read double from string '"+str+"'", e)
  }
}

/***********************************************/

Int String::toInt(const std::string &str)
{
  try
  {
    if(std::all_of(str.begin(), str.end(), ::isspace))
      return 0;

    return std::stoi(str);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("cannot read integer from string '"+str+"'", e)
  }
}

/***********************************************/

Bool String::contains(const std::string &str, const std::string &test)
{
  return str.size() && (str.find(test) != std::string::npos);
}

/***********************************************/

Bool String::startsWith(const std::string &str, const std::string &test)
{
  return (str.rfind(test, 0) == 0);
}

/***********************************************/

Bool String::endsWith(const std::string &str, const std::string &test)
{
  return (str.size() >= test.size()) && (str.compare(str.size()-test.size(), test.size(), test) == 0);
}

/***********************************************/

std::vector<std::string> String::split(const std::string &str, Char seperator)
{
  try
  {
    std::vector<std::string> parts;
    UInt pos1 = 0;
    UInt pos2 = 0;
    for(;;)
    {
      pos2 = str.find(seperator, pos2);
      parts.push_back(str.substr(pos1, pos2-pos1));
      if(pos2 == std::string::npos)
        break;
      pos1 = ++pos2;
    }
    return parts;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::vector<std::string> String::split(const std::string &str, const std::string &seperators)
{
  try
  {
    std::vector<std::string> parts;
    UInt pos1 = 0;
    UInt pos2 = 0;
    for(;;)
    {
      pos2 = str.find_first_of(seperators, pos2);
      parts.push_back(str.substr(pos1, pos2-pos1));
      if(pos2 == std::string::npos)
        break;
      pos1 = ++pos2;
    }
    return parts;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string String::replaceAll(const std::string &str, const std::string &search, const std::string &substitute)
{
  std::string result = str;
  UInt pos = 0;
  for(;;)
  {
    pos = result.find(search, pos);
    if(pos == std::string::npos)
        break;
    result.replace(pos, search.size(), substitute);
    pos += substitute.size();
  }
  return result;
}

/***********************************************/

std::string String::replaceAll(const std::string &str, const std::vector<std::pair<std::string, std::string>> &substitutes)
{
  std::string result = str;
  for(const auto &substitute : substitutes)
  {
    auto pos = result.find(substitute.first);
    while(pos != std::string::npos)
    {
      result.replace(pos, substitute.first.size(), substitute.second);
      pos = result.find(substitute.first, pos + substitute.second.size());
    }
  }
  return result;
}

/***********************************************/

std::regex String::wildcard2regex(const std::string &pattern)
{
  try
  {
    // escape ., +, \*, \? in regex pattern and
    // replace valid wildcards *, ? with regex equivalents .*, .
    std::string regexPattern = pattern;
    regexPattern = replaceAll(regexPattern, {{"\\*", "<backslash_star>"}, {"\\?", "<backslash_question_mark>"}});  // temp replacement
    regexPattern = replaceAll(regexPattern, {{"\\", "\\\\"}, {"^", "\\^"}, {".", "\\."}, {"$", "\\$"}, {"|", "\\|"}, {"(", "\\("},
                                             {")", "\\)"},   {"[", "\\["}, {"]", "\\]"}, {"+", "\\+"}, {"/", "\\/"}});
    regexPattern = replaceAll(regexPattern, {{"*", ".*"}, {"?", "."}});
    regexPattern = replaceAll(regexPattern, {{"<backslash_star>", "\\*"}, {"<backslash_question_mark>", "\\?"}}); // undo temp replacement
    return std::regex(regexPattern);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
