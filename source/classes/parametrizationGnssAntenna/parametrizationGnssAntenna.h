/***********************************************/
/**
* @file parametrizationGnssAntenna.h
*
* @brief Parametrization of antenna center variations.
*
* @author Torsten Mayer-Guerr
* @date 2012-11-08
*
*/
/***********************************************/

#ifndef __GROOPS_PARAMETRIZATIONGNSSANTENNA__
#define __GROOPS_PARAMETRIZATIONGNSSANTENNA__

// Latex documentation
#ifdef DOCSTRING_ParametrizationGnssAntenna
static const char *docstringParametrizationGnssAntenna = R"(
\section{ParametrizationGnssAntenna}\label{parametrizationGnssAntennaType}
Parametrization of antenna center variations. It will be used to set up the design matrix in a least squares adjustment.
Usually the paramtrization is setup separately for different \configClass{gnssType}{gnssType}.

If multiple parametrizations are given the parameters are sequently appended in the design matrix and parameter vector.
)";
#endif

/***********************************************/

#include "config/config.h"
#include "base/gnssType.h"
#include "gnss/gnss.h"

/**
* @defgroup parametrizationGnssAntennaGroup ParametrizationGnssAntenna
* @brief Parametrization of antenna center variations.
* @ingroup classesGroup
* The interface is given by @ref ParametrizationGnssAntenna.
* An Instance can be created by @ref readConfig. */
/// @{

/***** TYPES ***********************************/

class ParametrizationGnssAntenna;
class ParametrizationGnssAntennaBase;
typedef std::shared_ptr<ParametrizationGnssAntenna> ParametrizationGnssAntennaPtr;

/***** CLASS ***********************************/

/** @brief Parametrization of antenna center variations.
* An Instance of this class can be created by @ref readConfig. */
class ParametrizationGnssAntenna
{
  UInt parameterCount_;
  std::vector<UInt> index;
  std::vector<ParametrizationGnssAntennaBase*> base;

public:
  /// Constructor.
  ParametrizationGnssAntenna(Config &config, const std::string &name);

  /// Destructor.
  ~ParametrizationGnssAntenna();

  /** @brief Number of parameters.
  * This is the column count of the design matrix @a A. */
  UInt parameterCount() const {return parameterCount_;}

  /** @brief Name of parameters.
  * The names are appended to @a name. */
  void parameterName(std::vector<ParameterName> &name) const;

  /** @brief Partial Derivations of antenna center variations.
  * @return (1 x parameterCount) matrix with partial derivatives. */
  Matrix designMatrix(Angle azimut, Angle elevation);

  /** @brief creates an derived instance of this class. */
  static ParametrizationGnssAntennaPtr create(Config &config, const std::string &name) {return ParametrizationGnssAntennaPtr(new ParametrizationGnssAntenna(config, name));}
};

/***** FUNCTIONS *******************************/

/** @brief creates an instance of the class ParametrizationGnssAntenna.
* search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE.
* @param config the config node which includes the node with the options for this class
* @param name tag name in the config.
* @param var output: created class.
* @param mustSet if @a name is not found and @a mustSet=Config::MUSTSET, this function throws an exception instead of returning with FALSE.
* @param defaultValue ignored at the moment.
* @param annotation description of the function of this class.
* @relates ParametrizationGnssAntenna */
template<> Bool readConfig(Config &config, const std::string &name, ParametrizationGnssAntennaPtr &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***** CLASS ***********************************/

// Internal class
class ParametrizationGnssAntennaBase
{
public:
  virtual ~ParametrizationGnssAntennaBase() {}
  virtual UInt parameterCount() const = 0;
  virtual void parameterName(std::vector<ParameterName> &name) const = 0;
  virtual void designMatrix(Angle azimut, Angle elevation, MatrixSliceRef A) = 0;
};

/***********************************************/

#endif /* __GROOPS___ */
