/***********************************************/
/**
* @file loop.h
*
* @brief Loop.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2017-01-27
*
*/
/***********************************************/

#ifndef __GROOPS_LOOP__
#define __GROOPS_LOOP__

// Latex documentation
#ifdef DOCSTRING_Loop
static const char *docstringLoop = R"(
\section{Loop}\label{loopType}
Generates a sequence with variables to loop over.
The variable names can be set with \config{variableLoop...} and
the current values are assigned to the variables for each loop step.
See \reference{Loop and conditions}{general.loopsAndConditions} for usage.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"

/**
* @defgroup loopGroup Loop
* @brief Loop config elements of programs.
* @ingroup classesGroup
* The interface is given by @ref Loop.
* An Instance can be created by @ref readConfig. */
/// @{

/***** TYPES ***********************************/

class Loop;
typedef std::shared_ptr<Loop> LoopPtr;

/***** CLASS ***********************************/

/** @brief Loop.
* An instance of this class can be created with @ref readConfig. */
class Loop
{
public:
  /// Destructor.
  virtual ~Loop() {}

  /** @brief Returns the approximate total number of iterations. */
  virtual UInt count() const = 0;

  /** @brief Sets values of loop variables in @p varList for current iteration.
  * @return valid iteration step? */
  virtual Bool iteration(VariableList &varList) = 0;

  /** @brief creates an derived instance of this class. */
  static LoopPtr create(Config &config, const std::string &name);
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class Loop.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a loop is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] loop Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates Loop */
template<> Bool readConfig(Config &config, const std::string &name, LoopPtr &loop, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

#endif /* __GROOPS_DERIVATION__ */
