/***********************************************/
/**
* @file stringParser.cpp
*
* @brief string manipulation parser
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2016-03-26
*
*/
/***********************************************/

#include "base/importStd.h"
#include "base/constants.h"
#include "stringParser.h"
#include "expressionParser.h"
#include <regex>

/***********************************************/

static std::string parseUntil(const std::string &text, const char *search, const VariableList &varList, std::string::size_type &pos, Bool &resolved);
static std::string parseVariable(const std::string &text, const VariableList &varList, std::string::size_type &pos, Bool &resolved);

/***********************************************/

std::string parseUntil(const std::string &text, const char *search, const VariableList &varList, std::string::size_type &pos, Bool &resolved)
{
  std::string textPart;
  for(;;)
  {
    std::string::size_type posOld = pos;
    pos = text.find_first_of(search, pos);
    textPart += text.substr(posOld, pos-posOld);
    if((pos != std::string::npos) && (pos > 1) && (text.at(pos-1) == '#')) // escape '#' before?
      textPart.back() = text.at(pos++);                                    // replace '#' by excapted character
    else if((pos != std::string::npos) && (text.at(pos) == '{'))
      textPart += parseVariable(text, varList, ++pos, resolved);
    else
      break;
  }
  return textPart;
}

/***********************************************/

// parse after '{'
static std::string parseVariable(const std::string &text, const VariableList &varList, std::string::size_type &pos, Bool &resolved)
{
  Bool resolvedVar = TRUE;
  std::string textPart = parseUntil(text, "{:/}", varList, pos, resolvedVar);
  if(pos == std::string::npos)
    throw(Exception("missing closing '}'"));
  auto c = text.at(pos++);

  // textPart is variable
  // --------------------
  if(c == '}')
  {
    auto variable = varList.find(textPart);
    if(variable && resolvedVar)
      return variable->getParsedText(varList, resolved);
    resolved = FALSE;
    return '{'+textPart+'}';
  }

  // textPart is {expression:format}
  // -------------------------------
  if(c == ':')
  {
    // format string
    std::string format = parseUntil(text, "{}", varList, pos, resolvedVar);
    if(pos++ == std::string::npos)
      throw(Exception("missing closing '}' of {expression:format}"));
    if(resolvedVar)
    {
      try {return ExpressionVariable::parse(textPart, varList) % format;} // parse expression -> caluclate new result string
      catch(std::exception &/*e*/) {}
    }
    resolved = FALSE;
    return "{"+textPart+":"+format+"}";
  }

  // regex expression {text/regex/replace}
  // -------------------------------------
  if(c == '/')
  {
    std::string regexText = parseUntil(text, "{/", varList, pos, resolvedVar);
    if(pos++ == std::string::npos)
      throw(Exception("missing second '/' of {text/regex/replace}"));
    std::string replace = parseUntil(text, "{}", varList, pos, resolvedVar);
    if(pos++ == std::string::npos)
      throw(Exception("missing closing '}' of {text/regex/replace}"));
    if(!resolvedVar)
    {
      resolved = FALSE;
      return '{'+textPart+'/'+regexText+'/'+replace+'}';
    }

    // Escape sequences
    // \l lowercase next char
    // \u uppercase next char
    // \L lowercase until \E
    // \U uppercase until \E
    // \Q quote (disable) pattern metacharacters until \E
    // \E end either case modification or quoted section
    std::regex regex(regexText);
    std::string result;
    std::string::size_type posReplace = 0;
    Char c = 'E'; // end: normal mode
    for(;;)
    {
      std::string::size_type posOld = posReplace;
      posReplace = replace.find('\\', posOld);

      if(c != 'Q')
      {
        std::string resultPart = std::regex_replace(textPart, regex, replace.substr(posOld, posReplace-posOld));
        if(c == 'L') std::for_each(resultPart.begin(), resultPart.end(), [](auto &c){c=std::tolower(c);});
        if(c == 'U') std::for_each(resultPart.begin(), resultPart.end(), [](auto &c){c=std::toupper(c);});
        if(c == 'l' && !resultPart.empty()) resultPart.front() = std::tolower(resultPart.front());
        if(c == 'u' && !resultPart.empty()) resultPart.front() = std::toupper(resultPart.front());
        result += resultPart;
      }
      else
        result += replace.substr(posOld, posReplace-posOld);

      if(posReplace++ == std::string::npos) return result;      // end of string
      auto cNew = (posReplace < replace.size()) ? replace.at(posReplace) : ' ';
      if(((c != 'Q') && ((cNew == 'L') || (cNew == 'U') || (cNew == 'l') || (cNew == 'u') || (cNew == 'Q'))) || // not in quote mode  -> accept Q,L,U,l,u
         ((c != 'E') && (cNew == 'E')))                                                                         // not in normal mode -> accept additionally E
        c = replace.at(posReplace++);
      else
        result += "\\";
    }

   return result;
  }

  throw(Exception("unknown error"));
}

/***********************************************/

std::string StringParser::parse(const std::string &name, const std::string &text, const VariableList &varList, Bool &resolved)
{
  std::string result;
  std::string::size_type pos = 0;
  try
  {
    return parseUntil(text, "{", varList, pos, resolved);
  }
  catch(std::exception &e)
  {
    throw(Exception("Parser error in ("+name+"='"+text+"'), column="+(pos%"%i"s)+":\n"+e.what()));
  }
}

/***********************************************/

std::string StringParser::parse(const std::string &text, const VariableList &varList)
{
  try
  {
    Bool resolved = TRUE;
    std::string result = parse("(unknown)", text, varList, resolved);
    if(!resolved)
      throw(Exception("unresolved variables"));
    return result;
  }
  catch(std::exception &e)
  {
    throw(Exception("Parser error in '"+text+"'\n"+e.what()));
  }
}

/***********************************************/
