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

/***********************************************/

// parse after '{'
static std::string parseVar(const std::string &text, const VariableList &varList, std::string::size_type &pos, Bool &resolved)
{
  std::string result;
  Bool resolvedVar = TRUE;
  for(;;)
  {
    auto posTmp = text.find_first_of("{:}", pos);
    if(posTmp == std::string::npos)
      throw(Exception("missing closing '}'"));
    result += text.substr(pos, posTmp-pos);
    pos = posTmp;
    if(text.at(pos)!='{')
      break;
    // variable in variable
    pos++; // after the '{'
    result += parseVar(text, varList, pos, resolvedVar);
  }

  // result is variable?
  if(text.at(pos)=='}')
  {
    pos++;
    // trim eat white space
    auto start = result.find_first_not_of(" \t");
    if(start == std::string::npos)
      return ""; // empty string
    auto end = result.find_last_not_of(" \t");
    result = result.substr(start, end-start+1);

    auto variable = varList.find(result);
    if((!resolvedVar) || (!variable))
    {
      resolved = FALSE;
      return '{'+result+'}';
    }
    return variable->getParsedText(varList, resolved);
  }

  // result is expression
  pos++; // after the ':'

  // format string
  std::string format;
  for(;;)
  {
    auto posTmp = text.find_first_of("{}", pos);
    if(posTmp == std::string::npos)
      throw(Exception("missing closing '}'"));
    format += text.substr(pos, posTmp-pos);
    pos = posTmp;
    if(text.at(pos)!='{') // variable in format string
      break;
    pos++; // after the '{'
    // variable in variable
    format += parseVar(text, varList, pos, resolvedVar);
  }
  pos++; // after }

  if(!resolvedVar)
  {
    resolved = FALSE;
    return "{"+result+":"+format+"}";
  }

  // parse expression
  ExpressionPtr expr = Expression::parse(result);
  Double value = 0;
  try
  {
    value = expr->evaluate(varList);
  }
  catch(std::exception &/*e*/)
  {
    return "{"+result+":"+format+"}";
  }

  // parse format string -> caluclate new result string
  return value%format;
}

/***********************************************/

std::string StringParser::parse(const std::string &name, const std::string &text, const VariableList &varList, Bool &resolved)
{
  std::string::size_type pos = 0;

  try
  {
    resolved = TRUE;
    std::string result;
    for(;;)
    {
      std::string::size_type posOld = pos;
      pos = text.find('{', posOld);
      result += text.substr(posOld, pos-posOld);
      if(pos == std::string::npos)
        break;
      pos++;
      result += parseVar(text, varList, pos, resolved);
    }

    return result;
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
    Bool resolved;
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
