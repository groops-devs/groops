/***********************************************/
/**
* @file parametrizationTemporal.h
*
* @brief Parametrization of temporal variations.
*
* @author Torsten Mayer-Guerr
* @date 2011-12-28
*
*/
/***********************************************/

#ifndef __GROOPS_PARAMETRIZATIONTEMPORAL__
#define __GROOPS_PARAMETRIZATIONTEMPORAL__

// Latex documentation
#ifdef DOCSTRING_ParametrizationTemporal
static const char *docstringParametrizationTemporal = R"(
\section{ParametrizationTemporal}\label{parametrizationTemporalType}
This class gives a parametrization of time depending parameters (gravity field, positions, ...).
It will be used to set up the design matrix in a least squares adjustment.
If multiple parametrizations are given the coefficients in the parameter vector
are sequently appended.

Useally time intervals are defined half open meaning the last time belongs not to the interval.
This behaviour can be changed for the last interval with \config{includeLastTime}.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "base/parameterName.h"
#include "config/config.h"

/**
* @defgroup parametrizationTemporalGroup ParametrizationTemporal
* @brief Parametrization of temporal variations.
* @ingroup classesGroup
* The interface is given by @ref ParametrizationTemporal. */
/// @{

/***** TYPES ***********************************/

class ParametrizationTemporal;
class ParametrizationTemporalBase;
typedef std::shared_ptr<ParametrizationTemporal> ParametrizationTemporalPtr;

/***** CLASS ***********************************/

/** @brief Parametrization of temporal variations.
* An Instance of this class can be created by @ref readConfig. */
class ParametrizationTemporal
{
  UInt parameterCount_;
  std::vector<UInt> index;
  std::vector<ParametrizationTemporalBase*> representation;

public:
  /// Constructor.
  ParametrizationTemporal(Config &config, const std::string &name);

  /// Destructor.
  ~ParametrizationTemporal();

  /** @brief Estimate parameter in the given interval only.
  * Interval is defined as [@a timeStart, @a timeEnd).
  * Change result of @a parameterCount(), @a parameterName().
  * @return TRUE if parameters are changed */
  Bool setInterval(const Time &timeStart, const Time &timeEnd, Bool estimatePerArc=FALSE);

  /** @brief number of parameters. */
  UInt parameterCount() const {return parameterCount_;}

  /** @brief Multipliers of unknown parameters for a point in @a time. */
  Vector factors(const Time &time) const;

  /** @brief Multipliers of unknown parameters for a point in @a time.
  * @param time Time.
  * @param[out] index Position of multipliers in the design matrix.
  * @param[out] factor Multipliers of unknown parameters. */
  void factors(const Time &time, std::vector<UInt> &index, std::vector<Double> &factor) const;

  /** @brief Temporal design matrix.
  * The temporal design matrix @a A is constructed
  * from the basic design matrix @a B at this point in time. */
  void designMatrix(const Time &time, const const_MatrixSlice &B, MatrixSliceRef A) const;

  /** @brief Name of parameters.
  * The names are appended to @a name. */
  void parameterName(std::vector<ParameterName> &name) const;

  /** @brief Name of parameters.
  * @a baseName contain the name of each temporal component (e.g. x,y,z).
  * They are connected with ':', e.g. x:cos, y:cos, z:cos, x:sin, y:sin, z:sin.
  * The names are appended to @a name. */
  void parameterName(const std::vector<ParameterName> &baseName, std::vector<ParameterName> &name) const;

  /** @brief creates an derived instance of this class. */
  static ParametrizationTemporalPtr create(Config &config, const std::string &name) {return ParametrizationTemporalPtr(new ParametrizationTemporal(config, name));}
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class ParametrizationTemporal.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a representation without parameters is created.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] representation Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates ParametrizationTemporal */
template<> Bool readConfig(Config &config, const std::string &name, ParametrizationTemporalPtr &representation, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***** CLASS ***********************************/

// Internal class
class ParametrizationTemporalBase
{
public:
virtual ~ParametrizationTemporalBase() {}
virtual Bool setInterval(const Time &timeStart, const Time &timeEnd, Bool estimatePerArc) = 0;
virtual UInt parameterCount() const = 0;
virtual void factors(const Time &time, UInt startIndex, std::vector<UInt> &index, std::vector<Double> &factor) const = 0;
virtual void parameterName(std::vector<ParameterName> &name) const = 0;
};

/***********************************************/

#endif
