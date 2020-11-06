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
  TimeSeriesPtr     timeSeriesPtr;

public:
  LoopTimeSeries(Config &config);

  UInt count() const override;
  void init(VariableList &varList) override;
  void setValues(VariableList &varList) override;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline LoopTimeSeries::LoopTimeSeries(Config &config) : Loop()
{
  try
  {
    readConfig(config, "timeSeries",        timeSeriesPtr, Config::MUSTSET,   "",         "loop is called for every point in time");
    readConfig(config, "variableLoopTime",  nameTime,      Config::OPTIONAL,  "loopTime", "variable with time of each loop");
    readConfig(config, "variableLoopIndex", nameIndex,     Config::OPTIONAL,  "",         "variable with index of current iteration (starts with zero)");
    readConfig(config, "variableLoopCount", nameCount,     Config::OPTIONAL,  "",         "variable with total number of iterations");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline UInt LoopTimeSeries::count() const
{
  return times.size();
}

/***********************************************/

inline void LoopTimeSeries::init(VariableList &varList)
{
  times = timeSeriesPtr->times();

  if(index == NULLINDEX)
  {
    if(!nameTime.empty())  addVariable(nameTime,  varList);
    if(!nameIndex.empty()) addVariable(nameIndex, varList);
    if(!nameCount.empty()) addVariable(nameCount, varList);
  }
  index = 0;
}

/***********************************************/

inline void LoopTimeSeries::setValues(VariableList &varList)
{
  if(!nameTime.empty())  varList[nameTime]->setValue(times.at(index).mjd());
  if(!nameIndex.empty()) varList[nameIndex]->setValue(index);
  if(!nameCount.empty()) varList[nameCount]->setValue(count());
}

/***********************************************/

#endif
