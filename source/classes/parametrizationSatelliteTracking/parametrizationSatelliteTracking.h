/***********************************************/
/**
* @file parametrizationSatelliteTracking.h
*
* @brief SST parameter.
*
* @author Beate Klinger
* @date 2015-04-27
*
*/
/***********************************************/

#ifndef __GROOPS_PARAMETRIZATIONSATELLITETRACKING__
#define __GROOPS_PARAMETRIZATIONSATELLITETRACKING__

// Latex documentation
#ifdef DOCSTRING_ParametrizationSatelliteTracking
static const char *docstringParametrizationSatelliteTracking = R"(
\section{ParametrizationSatelliteTracking}\label{parametrizationSatelliteTrackingType}
This class defines parameters of Satellite-to-Satellite tracking observations.
It will be used to set up the design matrix in a least squares adjustment.
If multiple parametrizations are given the coefficients in the parameter vector
are sequently appended.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "base/parameterName.h"
#include "config/config.h"

/**
* @defgroup parametrizationSatelliteTrackingGroup ParametrizationSatelliteTracking
* @brief SST parameters.
* @ingroup classesGroup
* The interface is given by @ref ParametrizationSatelliteTracking. */
/// @{

/***** TYPES ***********************************/

class ParametrizationSatelliteTracking;
class ParametrizationSatelliteTrackingBase;
typedef std::shared_ptr<ParametrizationSatelliteTracking> ParametrizationSatelliteTrackingPtr;

/***** CLASS ***********************************/

/** @brief SST parameters.
* Parametrization of the (time variable) gravity field.
* Creates the observation equations of SST parameters (design matrix).
* An Instance of this class can be created by @ref readConfig. */
class ParametrizationSatelliteTracking
{
  UInt parameterCountA, parameterCountB;
  std::vector<UInt> indexA, indexB;
  std::vector<ParametrizationSatelliteTrackingBase*> parameter;

  void computeIndicies();

public:
  /// Constructor.
  ParametrizationSatelliteTracking(Config &config, const std::string &name);

  /// Destructor.
  ~ParametrizationSatelliteTracking();

  /** @brief Estimate parameter in the given interval only.
  * Change result of @a parameterCount(), @a parameterName().
  * @return TRUE if parameters are changed */
  Bool setInterval(const Time &timeStart, const Time &timeEnd);

  /** @brief Initialize arc.
  * Must be called at the beginning of every new arc.
  * @return TRUE if parameters are changed */
  Bool setIntervalArc(const Time &timeStart, const Time &timeEnd);

  /** @brief Number of unknown parameters.
  * This is the column count of the design matrix @a A. */
  UInt parameterCount() const {return parameterCountA;}

  /** @brief Number of unknown arc related parameters.
  * This is the column count of the design matrix @a B. */
  UInt parameterCountArc() const {return parameterCountB;}

  /** @brief Name of parameters.
  * The names are appended to @a name. */
  void parameterName(std::vector<ParameterName> &name) const;

  /** @brief Name of arc related parameters.
  * The names are appended to @a name. */
  void parameterNameArc(std::vector<ParameterName> &name) const;

  /** @brief Partial Derivations of force function.
  * @param sstType   SST observation type (0 = range, 1 = range rate, 2 = range acceleration)
  * @param time      GPS time
  * @param sst0      reference SST observations (Taylor point)
  * @param position1 (3*epochCount) in CRF [m]
  * @param position2 (3*epochCount) in CRF [m]
  * @param velocity1 (3*epochCount) in CRF [m/s]
  * @param velocity2 (3*epochCount) in CRF [m/s]
  * @param rotSat1   Sat1 -> CRF
  * @param rotSat2   Sat2 -> CRF
  * @param[out] A    (epochCount x parameterCount) matrix with partial derivatives
  * @param[out] B    (epochCount x parameterCountArc) matrix with partial derivatives */
  void compute(UInt sstType, const std::vector<Time> &time, const Vector &sst0,
              const Vector &position1, const Vector &position2, const Vector &velocity1, const Vector &velocity2,
              const std::vector<Rotary3d> &rotSat1, const std::vector<Rotary3d> &rotSat2,
              MatrixSliceRef A, MatrixSliceRef B);

  /** @brief creates an derived instance of this class. */
  static ParametrizationSatelliteTrackingPtr create(Config &config, const std::string &name) {return ParametrizationSatelliteTrackingPtr(new ParametrizationSatelliteTracking(config, name));}
};


/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class ParametrizationSatelliteTracking.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a parametrizationSatelliteTracking without parameters is created.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] parametrizationSatelliteTracking Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates ParametrizationSatelliteTracking */
template<> Bool readConfig(Config &config, const std::string &name, ParametrizationSatelliteTrackingPtr &parametrizationSatelliteTracking, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***** CLASS ***********************************/

// Internal class
class ParametrizationSatelliteTrackingBase
{
public:
virtual ~ParametrizationSatelliteTrackingBase() {}
virtual Bool isPerArc() const = 0;
virtual Bool setInterval(const Time &timeStart, const Time &timeEnd) = 0;
virtual UInt parameterCount()    const = 0;
virtual void parameterName(std::vector<ParameterName> &name) const = 0;
virtual void compute(UInt sstType, const std::vector<Time> &time, const Vector &sst0,
                     const Vector &position1, const Vector &position2, const Vector &velocity1, const Vector &velocity2,
                     const std::vector<Rotary3d> &rotSat1, const std::vector<Rotary3d> &rotSat2, MatrixSliceRef A) = 0;
};

/***********************************************/

#endif
