/***********************************************/
/**
* @file parametrizationAcceleration.h
*
* @brief Orbit force parameters.
*
* @author Torsten Mayer-Guerr
* @date 2014-03-18
*
*/
/***********************************************/

#ifndef __GROOPS_PARAMETRIZATIONACCELERATION__
#define __GROOPS_PARAMETRIZATIONACCELERATION__

// Latex documentation
#ifdef DOCSTRING_ParametrizationAcceleration
static const char *docstringParametrizationAcceleration = R"(
\section{ParametrizationAcceleration}\label{parametrizationAccelerationType}
This class defines parameters of satellite accelerations.
It will be used to set up the design matrix in a least squares adjustment.
If multiple parametrizations are given the coefficients in the parameter vector
are sequently appended.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "base/parameterName.h"
#include "config/config.h"
#include "files/fileSatelliteModel.h"
#include "classes/ephemerides/ephemerides.h"

/**
* @defgroup parametrizationAcceleration ParametrizationAcceleration
* @brief Orbit force parameters.
* @ingroup classesGroup
* The interface is given by @ref ParametrizationAcceleration.
* An Instance can be created by @ref readConfig. */
/// @{

/***** TYPES ***********************************/

class ParametrizationAcceleration;
class ParametrizationAccelerationBase;
typedef std::shared_ptr<ParametrizationAcceleration> ParametrizationAccelerationPtr;

/***** CLASS ***********************************/

/** @brief Orbit force parameters.
* Creates the observation equations of orbit parameters (design matrix).
* An Instance of this class can be created by @ref readConfig. */
class ParametrizationAcceleration
{
  UInt parameterCountA, parameterCountB;
  std::vector<UInt> indexA, indexB;
  std::vector<ParametrizationAccelerationBase*> parameter;

  void computeIndicies();

public:
  /// Constructor.
  ParametrizationAcceleration(Config &config, const std::string &name);

  /// Destructor.
 ~ParametrizationAcceleration();

  /** @brief Estimate parameter in the given interval only.
  * Change result of @a parameterCount(), @a parameterName().
  * @return TRUE if parameters are changed */
  Bool setInterval(const Time &timeStart, const Time &timeEnd);

  /** @brief Set the interval for the estimation arc related parameters.
  * Must be called at the beginning of every new arc.
  * Can change result of @a parameterCountArc(), @a parameterNameArc().
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
  * @param satellite Macro model.
  * @param time GPS time
  * @param position in CRF [m]
  * @param velocity in CRF [m/s]
  * @param rotSat   Sat -> CRF
  * @param rotEarth CRF -> TRF
  * @param ephemerides Position of Sun and Moon
  * @param[out] A (3 x parameterCount) matrix with partial derivatives in TRF(!) [m/s^2]
  * @param[out] B (3 x parameterCountArc) matrix with partial derivatives in TRF(!) [m/s^2] */
  void compute(SatelliteModelPtr satellite, const Time &time, const Vector3d &position, const Vector3d &velocity,
               const Rotary3d &rotSat, const Rotary3d &rotEarth, EphemeridesPtr ephemerides, MatrixSliceRef A, MatrixSliceRef B);

  /** @brief creates an derived instance of this class. */
  static ParametrizationAccelerationPtr create(Config &config, const std::string &name) {return ParametrizationAccelerationPtr(new ParametrizationAcceleration(config, name));}
};


/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class ParametrizationAcceleration.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a parametrizationAcceleration without parameters is created.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] parametrizationAcceleration Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates ParametrizationAcceleration */
template<> Bool readConfig(Config &config, const std::string &name, ParametrizationAccelerationPtr &parametrizationAcceleration, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***** CLASS ***********************************/

// Internal class
class ParametrizationAccelerationBase
{
public:
virtual ~ParametrizationAccelerationBase() {}
virtual Bool isPerArc() const = 0;
virtual Bool setInterval(const Time &timeStart, const Time &timeEnd) = 0;
virtual UInt parameterCount() const = 0;
virtual void parameterName(std::vector<ParameterName> &name) const = 0;
virtual void compute(SatelliteModelPtr satellite, const Time &time, const Vector3d &position, const Vector3d &velocity,
                     const Rotary3d &rotSat, const Rotary3d &rotEarth, EphemeridesPtr ephemerides, MatrixSliceRef A) = 0;
};

/***********************************************/

#endif
