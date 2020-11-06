/***********************************************/
/**
* @file parameterName.cpp
*
* @brief Parameter name representation.
*
* @author Sebastian Strasser
* @date 2017-05-23
*
*/
/***********************************************/

#include "base/importStd.h"
#include "base/string.h"
#include "base/parameterName.h"

/***********************************************/

ParameterName::ParameterName(const std::string &object, const std::string &type, const std::string &temporal, const Time &timeStart, const Time &timeEnd)
  : object(object), type(type), temporal(temporal)
{
  try
  {
    if(timeStart != Time())
    {
      std::stringstream ss;
      ss<<timeStart.dateTimeStr();
      if(timeEnd != Time())
        ss<<"_"<<timeEnd.dateTimeStr();
      interval = ss.str();
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

ParameterName::ParameterName(const std::string &object, const std::string &type, const std::string &temporal, const std::string &interval)
  : object(object), type(type), temporal(temporal), interval(interval)
{
}

/***********************************************/

ParameterName ParameterName::fromStr(const std::string &str)
{
  try
  {
    std::vector<std::string> parts = String::split(str, sep);

    ParameterName parameterName;
    if(parts.size() > 0) parameterName.object   = parts.at(0);
    if(parts.size() > 1) parameterName.type     = parts.at(1);
    if(parts.size() > 2) parameterName.temporal = parts.at(2);
    if(parts.size() > 3) parameterName.interval = parts.at(3);
    if(parts.size() > 4) throw(Exception("Parameter name contains more than four parts: "+str));

    return parameterName;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool ParameterName::combine(const ParameterName &name)
{
  try
  {
    if(!fuzzyMatch(name))
      return FALSE;

    if(object.empty())   object   = name.object;
    if(type.empty())     type     = name.type;
    if(interval.empty()) interval = name.interval;
    if(temporal.empty()) temporal = name.temporal;

    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool ParameterName::fuzzyMatch(const ParameterName &other) const
{
  try
  {
    return (object   == other.object   || object.empty()   || other.object.empty())   &&
           (type     == other.type     || type.empty()     || other.type.empty())     &&
           (interval == other.interval || interval.empty() || other.interval.empty()) &&
           (temporal == other.temporal || temporal.empty() || other.temporal.empty());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool ParameterName::operator==(const ParameterName &other) const
{
  try
  {
    return (object == other.object) && (type == other.type ) && (interval == other.interval ) && (temporal == other.temporal);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool ParameterName::operator<(const ParameterName &other) const
{
  try
  {
    Int compareObject = object.compare(other.object);
    if(compareObject < 0) return TRUE;
    if(compareObject > 0) return FALSE;

    Int compareType = type.compare(other.type);
    if(compareType < 0) return TRUE;
    if(compareType > 0) return FALSE;

    Int compareInterval = interval.compare(other.interval);
    if(compareInterval < 0) return TRUE;
    if(compareInterval > 0) return FALSE;

    Int compareTemporal = temporal.compare(other.temporal);
    if(compareTemporal < 0) return TRUE;

    return FALSE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
