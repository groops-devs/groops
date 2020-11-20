/***********************************************/
/**
* @file loopLoop.h
*
* @brief Loop over nested loops.
*
* @author Matthias Ellmer
* @author Sebastian Strasser
* @date 2017-07-10
*
*/
/***********************************************/

#ifndef __GROOPS_LOOPLOOP__
#define __GROOPS_LOOPLOOP__

// Latex documentation
#ifdef DOCSTRING_Loop
static const char *docstringLoopLoop = R"(
\subsection{Loop}
Loop over nested loops. First \config{loop} is outermost loop, every subsequent \config{loop} is one level below the previous \config{loop}.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "classes/loop/loop.h"

/***** CLASS ***********************************/

/** @brief Loop over nested loops.
* @ingroup LoopGroup
* @see Loop */
class LoopLoop : public Loop
{
  std::vector<LoopConfig> loopConfigs;
  std::vector<LoopPtr>    loops;
  UInt                    index;
  std::string             nameIndex;

public:
  LoopLoop(Config &config);

  UInt count() const override;
  Bool iteration(VariableList &varList) override;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline LoopLoop::LoopLoop(Config &config)
{
  try
  {
    readConfig(config, "loop",               loopConfigs, Config::MUSTSET,  "", "subloop");
    readConfig(config, "variableLoopIndex",  nameIndex,   Config::OPTIONAL, "", "variable with index of current iteration (starts with zero)");
    if(isCreateSchema(config))
      return;

    loops.resize(loopConfigs.size());
    index = 0;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline UInt LoopLoop::count() const
{
  UInt count = 1;
  for(auto loop : loops)
    if(loop)
      count *= loop->count();
  return std::max(index, count);
}

/***********************************************/

inline Bool LoopLoop::iteration(VariableList &varList)
{
  try
  {
    std::function<Bool(UInt)> initLoop = [&](UInt i) -> Bool
    {
      loops.at(i) = loopConfigs.at(i).read(varList);
      while(loops.at(i)->iteration(varList))
        if((i+1 >= loops.size()) || initLoop(i+1))
          return TRUE;
      return FALSE;
    };

    if(index == 0)
    {
      if(!nameIndex.empty())
        addVariable(nameIndex, index, varList);
      index++;
      return initLoop(0);
    }

    for(UInt i=loops.size(); i-->0;)
      while(loops.at(i)->iteration(varList))
        if((i+1 >= loops.size()) || initLoop(i+1))
        {
          if(!nameIndex.empty())
            varList[nameIndex]->setValue(index);
          index++;
          return TRUE;
        }

    return FALSE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
