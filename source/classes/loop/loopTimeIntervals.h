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

#ifndef __GROOPS_LOOPTIMEINTERVALS__
#define __GROOPS_LOOPTIMEINTERVALS__

// Latex documentation
#ifdef DOCSTRING_Loop
static const char *docstringLoopTimeIntervals = R"(
\subsection{TimeIntervals}\label{loopType:timeIntervals}
Loop over the intervals between points in time.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "classes/loop/loop.h"
#include "classes/timeSeries/timeSeries.h"

/***** CLASS ***********************************/

/** @brief Loop over time intervals.
* @ingroup loopGroup
* @see Loop */
class LoopTimeIntervals : public Loop
{
  std::string       nameTimeStart, nameTimeEnd, nameIndex, nameCount;
  std::vector<Time> times;

public:
  LoopTimeIntervals(Config &config);

  UInt count() const override {return times.size()-1;}
  Bool iteration(VariableList &varList) override;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline LoopTimeIntervals::LoopTimeIntervals(Config &config)
{
  try
  {
    TimeSeriesPtr timesIntervalPtr;

    readConfig(config, "timeIntervals",         timesIntervalPtr, Config::MUSTSET,   "",            "loop is called for every interval");
    readConfig(config, "variableLoopTimeStart", nameTimeStart,    Config::OPTIONAL,  "loopTime",    "variable with starting time of each interval");
    readConfig(config, "variableLoopTimeEnd",   nameTimeEnd,      Config::OPTIONAL,  "loopTimeEnd", "variable with ending time of each interval");
    readConfig(config, "variableLoopIndex",     nameIndex,        Config::OPTIONAL,  "",            "variable with index of current iteration (starts with zero)");
    readConfig(config, "variableLoopCount",     nameCount,        Config::OPTIONAL,  "",            "variable with total number of iterations");
    readConfigCondition(config);
    if(isCreateSchema(config)) return;

    times = timesIntervalPtr->times();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Bool LoopTimeIntervals::iteration(VariableList &varList)
{
  if(index() >= count())
    return FALSE;

  if(!nameTimeStart.empty()) varList.setVariable(nameTimeStart, times.at(index()).mjd());
  if(!nameTimeEnd.empty())   varList.setVariable(nameTimeEnd,   times.at(index()+1).mjd());
  if(!nameIndex.empty())     varList.setVariable(nameIndex,     index());
  if(!nameCount.empty())     varList.setVariable(nameCount,     count());

  return checkCondition(varList);
}

/***********************************************/

#endif
