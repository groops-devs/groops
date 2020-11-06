/***********************************************/
/**
* @file loopUniformSampling.h
*
* @brief Loop over sequence of numbers.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2017-01-27
*
*/
/***********************************************/

#ifndef __GROOPS_LOOPUNIFORMSAMPLING__
#define __GROOPS_LOOPUNIFORMSAMPLING__

// Latex documentation
#ifdef DOCSTRING_Loop
static const char *docstringLoopUniformSampling = R"(
\subsection{UniformSampling}
Loop over sequence of numbers.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "classes/loop/loop.h"

/***** CLASS ***********************************/

/** @brief Loop over sequence of numbers.
* @ingroup LoopGroup
* @see Loop */
class LoopUniformSampling : public Loop
{
  std::string nameNumber, nameIndex, nameCount;
  std::vector<Double> numbers;
  Double rangeStart, rangeEnd, sampling;

public:
  LoopUniformSampling(Config &config);

  UInt count() const override;
  void init(VariableList &varList) override;
  void setValues(VariableList &varList) override;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline LoopUniformSampling::LoopUniformSampling(Config &config) : Loop()
{
  try
  {
    readConfig(config, "rangeStart",         rangeStart, Config::MUSTSET,   "",           "start of range");
    readConfig(config, "rangeEnd",           rangeEnd,   Config::MUSTSET,   "",           "end of range (inclusive)");
    readConfig(config, "sampling",           sampling,   Config::MUSTSET,   "",           "sampling");
    readConfig(config, "variableLoopNumber", nameNumber, Config::OPTIONAL,  "loopNumber", "name of the variable to be replaced");
    readConfig(config, "variableLoopIndex",  nameIndex,  Config::OPTIONAL,  "",           "variable with index of current iteration (starts with zero)");
    readConfig(config, "variableLoopCount",  nameCount,  Config::OPTIONAL,  "",           "variable with total number of iterations");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline UInt LoopUniformSampling::count() const
{
  return numbers.size();
}

/***********************************************/

inline void LoopUniformSampling::init(VariableList &varList)
{
  numbers.clear();
  for(UInt i=0; rangeStart+i*sampling<=rangeEnd; i++)
    numbers.push_back(rangeStart + i * sampling);

  if(index == NULLINDEX)
  {
    if(!nameNumber.empty()) addVariable(nameNumber, varList);
    if(!nameIndex.empty())  addVariable(nameIndex,  varList);
    if(!nameCount.empty())  addVariable(nameCount,  varList);
  }
  index = 0;
}

/***********************************************/

inline void LoopUniformSampling::setValues(VariableList &varList)
{
  if(!nameNumber.empty()) varList[nameNumber]->setValue(numbers.at(index));
  if(!nameIndex.empty())  varList[nameIndex]->setValue(index);
  if(!nameCount.empty())  varList[nameCount]->setValue(count());
}

/***********************************************/

#endif
