/***********************************************/
/**
* @file loopTimeSeries.h
*
* @brief Loop over points in time.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2017-01-27
*
*/
/***********************************************/

#ifndef __GROOPS_LOOPTEMPORAL__
#define __GROOPS_LOOPTEMPORAL__

// Latex documentation
#ifdef DOCSTRING_Loop
static const char *docstringLoopTimeSeries = R"(
\subsection{TimeSeries}
Loop over points in time.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "classes/loop/loop.h"
#include "classes/timeSeries/timeSeries.h"

/***** CLASS ***********************************/

/** @brief Loop over points in time.
* @ingroup LoopGroup
* @see Loop */
class LoopTimeSeries : public Loop
{
  std::string       nameTime, nameIndex, nameCount;
  std::vector<Time> times;

public:
  LoopTimeSeries(Config &config);

  UInt count() const override {return times.size();}
  Bool iteration(VariableList &varList) override;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline LoopTimeSeries::LoopTimeSeries(Config &config)
{
  try
  {
    TimeSeriesPtr timeSeriesPtr;

    readConfig(config, "timeSeries",        timeSeriesPtr, Config::MUSTSET,   "",         "loop is called for every point in time");
    readConfig(config, "variableLoopTime",  nameTime,      Config::OPTIONAL,  "loopTime", "variable with time of each loop");
    readConfig(config, "variableLoopIndex", nameIndex,     Config::OPTIONAL,  "",         "variable with index of current iteration (starts with zero)");
    readConfig(config, "variableLoopCount", nameCount,     Config::OPTIONAL,  "",         "variable with total number of iterations");
    readConfigCondition(config);
    if(isCreateSchema(config)) return;

    times = timeSeriesPtr->times();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Bool LoopTimeSeries::iteration(VariableList &varList)
{
  if(index() >= count())
    return FALSE;

  if(!nameTime.empty())  varList.setVariable(nameTime,  times.at(index()).mjd());
  if(!nameIndex.empty()) varList.setVariable(nameIndex, index());
  if(!nameCount.empty()) varList.setVariable(nameCount, count());

  return checkCondition(varList);
}

/***********************************************/

#endif
