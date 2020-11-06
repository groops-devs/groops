/***********************************************/
/**
* @file parameterSelector.h
*
* @brief Index vector from selected parameters.
*
* @author Sebastian Strasser
* @date 2018-05-08
*
*/
/***********************************************/

#ifndef __GROOPS_PARAMETERSELECTOR__
#define __GROOPS_PARAMETERSELECTOR__

// Latex documentation
#ifdef DOCSTRING_ParameterSelector
static const char *docstringParameterSelector = R"(
\section{ParameterSelector}\label{parameterSelectorType}
This class provides an index vector from selected parameters,
which can be used e.g. to reorder a normal equation matrix.
The size of the index vector determines the size of the new matrix.
Entries are the indices of the selected parameters in the provided
parameter list or NULLINDEX for zero/new parameters.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "base/parameterName.h"
#include "config/config.h"

/**
* @defgroup parameterSelector ParameterSelector
* @brief Index vector from selected parameters.
* @ingroup classesGroup
* The interface is given by @ref ParameterSelector.
* An Instance can be created by @ref readConfig. */
/// @{

/***** TYPES ***********************************/

class ParameterSelector;
class ParameterSelectorBase;
typedef std::shared_ptr<ParameterSelector> ParameterSelectorPtr;

/***** CLASS ***********************************/

/** @brief Index vector from selected parameters.
* This class provides an index vector from selected parameters.
* An Instance of this class can be created by @ref readConfig. */
class ParameterSelector
{
  VariableList varList;
  std::vector<ParameterSelectorBase*> parameters;

public:
  /// Constructor.
  ParameterSelector(Config &config, const std::string &name);

  /// Destructor.
 ~ParameterSelector();

  /** @brief Returns the index vector. */
  std::vector<UInt> indexVector(const std::vector<ParameterName> &parameterNames);

  /** @brief Returns an index vector containing all indexes from 0 to @p referenceLength that are not in @p vector. */
  static std::vector<UInt> indexVectorComplement(std::vector<UInt> vector, UInt referenceLength);

  /** @brief creates an derived instance of this class. */
  static ParameterSelectorPtr create(Config &config, const std::string &name) {return ParameterSelectorPtr(new ParameterSelector(config, name));}
};


/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class ParameterSelector.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a ParameterSelector with zero-size matrix is created.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] parameterSelector Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates ParameterSelector */
template<> Bool readConfig(Config &config, const std::string &name, ParameterSelectorPtr &parameterSelector, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***** CLASS ***********************************/

// Internal class
class ParameterSelectorBase
{
public:
  virtual ~ParameterSelectorBase() {}
  virtual std::vector<UInt> indexVector(const std::vector<ParameterName> &parameterNames, VariableList varList) = 0;
};

/***********************************************/

#endif
