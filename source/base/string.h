/***********************************************/
/**
* @file string.h
*
* @brief Miscellaneous string functions.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2019-09-04
*
*/
/***********************************************/

#ifndef __GROOPS_STRING__
#define __GROOPS_STRING__

/***********************************************/

#include "base/importStd.h"
#include <regex>

/***** CLASS ***********************************/

/**
* @brief Miscellaneous string functions.
* @ingroup base */
namespace String
{
  /** @brief The @p str is converted to upper case. */
  std::string upperCase(const std::string &str);

  /** @brief The @p str is converted to lower case. */
  std::string lowerCase(const std::string &str);

  /** @brief Remove all leading and trailing spaces from @p str. */
  std::string trim(const std::string &str);

  /** @brief Convert to Double. Returns 0 if substring is all white spaces. */
  Double toDouble(const std::string &str);

  /** @brief Convert to Int. Returns 0 if substring is all white spaces. */
  Int toInt(const std::string &str);

  /** @brief test wether the @a str starts with @p test. */
  Bool startsWith(const std::string &str, const std::string &test);

  /** @brief Replace all occurrences of @p search in @p str with the @p substitute string. */
  std::string replaceAll(const std::string &str, const std::string &search, const std::string &substitute);

  /** @brief Replace all occurrences of @p substitutes.first with the @p substitutes.second string. */
  std::string replaceAll(const std::string &str, const std::vector<std::pair<std::string, std::string>> &substitutes);

  /** @brief Returns a splitted string separated by @p seperator. */
  std::vector<std::string> split(const std::string &str, Char seperator);

  /** @brief Returns a splitted string separated by any character of @p seperators. */
  std::vector<std::string> split(const std::string &str, const std::string &seperators);

  /** @brief Convert a simple wildcard pattern into an equivalent regex expression.
  * Wildcards *, ? can be escaped \*, \?. */
  std::regex wildcard2regex(const std::string &pattern);
}

/***********************************************/

#endif /* __GROOPS___ */
