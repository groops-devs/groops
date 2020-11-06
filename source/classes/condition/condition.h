/***********************************************/
/**
* @file condition.h
*
* @brief Condition.
*
* @author Torsten Mayer-Guerr
* @date 2018-05-18
*
*/
/***********************************************/

#ifndef __GROOPS_CONDITION__
#define __GROOPS_CONDITION__

// Latex documentation
#ifdef DOCSTRING_Condition
static const char *docstringCondition = R"(
\section{Condition}\label{conditionType}
Test for conditions. See \reference{Loop and conditions}{general.loopsAndConditions} for usage.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"

/**
* @defgroup conditionGroup Condition
* @brief Conditional application of config elements or programs.
* @ingroup classesGroup
* The interface is given by @ref Condition.
* An Instance can be created by @ref readConfig. */
/// @{

/***** TYPES ***********************************/

class Condition;
typedef std::shared_ptr<Condition> ConditionPtr;

/***** CLASS ***********************************/

/** @brief Condition.
* An instance of this class can be created with @ref readConfig. */
class Condition
{
public:
  /// Destructor.
  virtual ~Condition() {}

  virtual Bool condition(const VariableList &varList) const = 0;

  /** @brief creates an derived instance of this class. */
  static ConditionPtr create(Config &config, const std::string &name);
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class Condition.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a condition is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] condition Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates Condition */
template<> Bool readConfig(Config &config, const std::string &name, ConditionPtr &condition, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

#endif /* __GROOPS_DERIVATION__ */
