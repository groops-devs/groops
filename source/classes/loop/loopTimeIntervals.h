/***********************************************/
/**
* @file loopTimeIntervals.h
*
* @brief Loop over time intervals.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2017-01-27
*
*/
/***********************************************/

#ifndef __GROOPS_LoopTimeIntervals__
#define __GROOPS_LoopTimeIntervals__

// Latex documentation
#ifdef DOCSTRING_Loop
static const char *docstringLoopTimeIntervals = R"(
\subsection{TimeIntervals}
Loop over the intervals between points in time.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "classes/loop/loop.h"
#include "classes/timeSeries/timeSeries.h"

/***** CLASS ***********************************/

/** @brief Loop over time intervals.
* @ingroup LoopGroup
* @see Loop */
class LoopTimeIntervals : public Loop
{
  std::string       nameTimeStart, nameTimeEnd, nameIndex, nameCount;
  std::vector<Time> timesInterval;
  TimeSeriesPtr     timesIntervalPtr;

public:
  LoopTimeIntervals(Config &config);

  UInt count() const override;
  void init(VariableList &varList) override;
  void setValues(VariableList &varList) override;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline LoopTimeIntervals::LoopTimeIntervals(Config &config) : Loop()
{
  try
  {
    readConfig(config, "timeIntervals",         timesIntervalPtr, Config::MUSTSET,   "",            "loop is called for every interval");
    readConfig(config, "variableLoopTimeStart", nameTimeStart,    Config::OPTIONAL,  "loopTime",    "variable with starting time of each interval");
    readConfig(config, "variableLoopTimeEnd",   nameTimeEnd,      Config::OPTIONAL,  "loopTimeEnd", "variable with ending time of each interval");
    readConfig(config, "variableLoopIndex",     nameIndex,        Config::OPTIONAL,  "",            "variable with index of current iteration (starts with zero)");
    readConfig(config, "variableLoopCount",     nameCount,        Config::OPTIONAL,  "",            "variable with total number of iterations");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline UInt LoopTimeIntervals::count() const
{
  return timesInterval.size()-1;
}

/***********************************************/

inline void LoopTimeIntervals::init(VariableList &varList)
{
  timesInterval = timesIntervalPtr->times();

  if(index == NULLINDEX)
  {
    if(!nameTimeStart.empty()) addVariable(nameTimeStart, varList);
    if(!nameTimeEnd.empty())   addVariable(nameTimeEnd,   varList);
    if(!nameIndex.empty())     addVariable(nameIndex,     varList);
    if(!nameCount.empty())     addVariable(nameCount,     varList);
  }
  index = 0;
}

/***********************************************/

inline void LoopTimeIntervals::setValues(VariableList &varList)
{
  if(!nameTimeStart.empty()) varList[nameTimeStart]->setValue(timesInterval.at(index).mjd());
  if(!nameTimeEnd.empty())   varList[nameTimeEnd]->setValue(timesInterval.at(index+1).mjd());
  if(!nameIndex.empty())     varList[nameIndex]->setValue(index);
  if(!nameCount.empty())     varList[nameCount]->setValue(count());
}

/***********************************************/

#endif
