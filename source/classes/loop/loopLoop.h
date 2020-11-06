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
   std::string nameIndex/*, nameCount*/;
   std::vector<LoopPtr> loops;

   void skipEmpty(VariableList varList);

public:
  LoopLoop(Config &config);

  UInt count() const override;
  void init(VariableList &varList) override;
  void setValues(VariableList &varList) override;
  void next(VariableList &varList) override;
  Bool finished() const override;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline LoopLoop::LoopLoop(Config &config) : Loop()
{
  try
  {
    readConfig(config, "loop", loops, Config::MUSTSET,   "", "subloop");
    readConfig(config, "variableLoopIndex",  nameIndex,  Config::OPTIONAL,  "", "variable with index of current iteration (starts with zero)");
    //readConfig(config, "variableLoopCount",  nameCount,  Config::OPTIONAL,  "", "variable with total number of iterations"); // should not be used since count is generated dynamically during iterations
    if(isCreateSchema(config))
      return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline UInt LoopLoop::count() const
{
  return index;
}

/***********************************************/

inline void LoopLoop::init(VariableList &varList)
{
  try
  {
    if(index == NULLINDEX)
    {
      if(!nameIndex.empty())  addVariable(nameIndex,  varList);
      //if(!nameCount.empty())  addVariable(nameCount,  varList);
    }
    index = 0;

    for(auto loop : loops)
    {
      loop->init(varList);
      if(loop->finished())
        break;
      loop->setValues(varList);
    }

    skipEmpty(varList);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void LoopLoop::setValues(VariableList &varList)
{
  for(auto loop : loops)
    loop->setValues(varList);

  if(!nameIndex.empty())  varList[nameIndex]->setValue(index);
  //if(!nameCount.empty())  varList[nameCount]->setValue(count());
}

/***********************************************/

inline void LoopLoop::next(VariableList &varList)
{
  try
  {
    loops.back()->next(varList);
    skipEmpty(varList);
    index++;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Bool LoopLoop::finished() const
{
  return loops.front()->finished();
}

/***********************************************/

inline void LoopLoop::skipEmpty(VariableList varList)
{
  try
  {
    while(!loops.front()->finished() && loops.back()->finished())
    {
      for(UInt i = loops.size(); i --> 0; )
      {
        loops.at(i)->next(varList);
        if(loops.at(i)->finished())
          continue;

        for(UInt j = i+1; j < loops.size(); j++)
        {
          loops.at(j-1)->setValues(varList);
          loops.at(j)->init(varList);
        }

        if(!loops.back()->finished())
          return;
        else
          break;
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
