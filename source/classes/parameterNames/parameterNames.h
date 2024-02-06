/***********************************************/
/**
* @file parameterNames.h
*
* @brief Generate parameter names.
*
* @author Torsten Mayer-Guerr
* @date 2020-05-29
*/
/***********************************************/

#ifndef __GROOPS_PARAMETERNAMES__
#define __GROOPS_PARAMETERNAMES__

// Latex documentation
#ifdef DOCSTRING_ParameterNames
static const char *docstringParameterNames = R"(
\section{ParameterNames}\label{parameterNamesType}
Generate a list of parameter names. All parameters are appended.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "base/parameterName.h"

/**
* @defgroup parameterNamesGroup Parameter names
* @brief Generate parameter names.
* @ingroup classesGroup
* The interface is given by @ref ParameterNames. */
/// @{

/***** TYPES ***********************************/

class ParameterNames;
typedef std::shared_ptr<ParameterNames> ParameterNamesPtr;

/***** CLASS ***********************************/

/** @brief Generate parameter names.
* An Instance of this class can be created by @ref readConfig. */
class ParameterNames
{
  std::vector<ParameterName> names;

public:
  /// Constructor.
  ParameterNames(Config &config, const std::string &name);

  /** @brief Point distribution. */
  const std::vector<ParameterName> &parameterNames() const {return names;}

  /** @brief creates an derived instance of this class. */
  static ParameterNamesPtr create(Config &config, const std::string &name) {return ParameterNamesPtr(new ParameterNames(config, name));}
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class ParameterNames.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and an class without names is created.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] var Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates ParameterNames */
template<> Bool readConfig(Config &config, const std::string &name, ParameterNamesPtr &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***** CLASS ***********************************/

// Internal class
class ParameterNamesBase
{
public:
  std::vector<ParameterName> names;

  virtual ~ParameterNamesBase() {}
};

/***********************************************/

#endif /* __GROOPS_GRID__ */
