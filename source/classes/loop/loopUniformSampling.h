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
  std::string         nameNumber, nameIndex, nameCount;
  std::vector<Double> numbers;
  UInt                index;

public:
  LoopUniformSampling(Config &config);

  UInt count() const override {return numbers.size();}
  Bool iteration(VariableList &varList) override;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline LoopUniformSampling::LoopUniformSampling(Config &config)
{
  try
  {
    Double rangeStart, rangeEnd, sampling;

    readConfig(config, "rangeStart",         rangeStart, Config::MUSTSET,   "",           "start of range");
    readConfig(config, "rangeEnd",           rangeEnd,   Config::MUSTSET,   "",           "end of range (inclusive)");
    readConfig(config, "sampling",           sampling,   Config::MUSTSET,   "",           "sampling");
    readConfig(config, "variableLoopNumber", nameNumber, Config::OPTIONAL,  "loopNumber", "name of the variable to be replaced");
    readConfig(config, "variableLoopIndex",  nameIndex,  Config::OPTIONAL,  "",           "variable with index of current iteration (starts with zero)");
    readConfig(config, "variableLoopCount",  nameCount,  Config::OPTIONAL,  "",           "variable with total number of iterations");
    if(isCreateSchema(config)) return;

    for(UInt i=0; rangeStart+i*sampling<=rangeEnd; i++)
      numbers.push_back(rangeStart + i*sampling);
    index = 0;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Bool LoopUniformSampling::iteration(VariableList &varList)
{
  if(index >= count())
    return FALSE;

  if(!nameNumber.empty()) varList.setVariable(nameNumber, numbers.at(index));
  if(!nameIndex.empty())  varList.setVariable(nameIndex,  index);
  if(!nameCount.empty())  varList.setVariable(nameCount,  count());

  index++;
  return TRUE;
}

/***********************************************/

#endif
