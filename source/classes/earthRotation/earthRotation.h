/***********************************************/
/**
* @file earthRotation.h
*
* @brief Transformation between CRF and TRF.
*
* @author Torsten Mayer-Guerr
* @date 2001-11-12
*
*/
/***********************************************/

#ifndef __GROOPS_EARTHROTATION__
#define __GROOPS_EARTHROTATION__

// Latex documentation
#ifdef DOCSTRING_EarthRotation
static const char *docstringEarthRotation = R"(
\section{EarthRotation}\label{earthRotationType}
This class realize the transformation between a terestrial
reference frame (TRF) and a celestial reference frame (CRF).
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"

/**
* @defgroup earthRotationGroup EarthRotation
* @brief Transformation between CRF and TRF.
* @ingroup classesGroup
* The interface is given by @ref EarthRotation. */
/// @{

/***** TYPES ***********************************/

class EarthRotation;
typedef std::shared_ptr<EarthRotation> EarthRotationPtr;

/***** CLASS ***********************************/

/** @brief Transformation between CRF and TRF.
* An Instance of this class can be created by @ref readConfig. */
class EarthRotation
{
public:
  /// Destructor.
  virtual ~EarthRotation() {}

  /** @brief Rotary matrix.
  * Inertial system (CRF) -> earth fixed system (TRF).
  * @param timeGPS modified julian date (MJD) in GPS time system. */
  virtual Rotary3d rotaryMatrix(const Time &timeGPS) const;

  /** @brief Instantaneous rotation vector of Earth rotation.
  * Contains the complete rotation with precession, nutation, polar wobble.
  * Given in inertial system (CRF). [rad/s].
  * @param timeGPS modified julian date (MJD) in GPS time system. */
  virtual Vector3d rotaryAxis(const Time &timeGPS) const;

  /** @brief Time derivative of the instantaneous rotation vector.
  * Given in inertial system (CRF). [rad/s^2].
  * @param timeGPS modified julian date (MJD) in GPS time system. */
  virtual Vector3d rotaryAxisDerivate(const Time &timeGPS)  const;

  /** @brief Earth orientation parameter.
  * All models and corrections applied.
  * @param timeGPS modified julian date (MJD) in GPS time system.
  * @param[out] xp polar motion [rad]
  * @param[out] yp polar motion [rad]
  * @param[out] sp TIO locator [rad]
  * @param[out] deltaUT UT1-UTC [seconds]
  * @param[out] LOD length of day [seconds]
  * @param[out] X precession & nutation [rad]
  * @param[out] Y precession & nutation [rad]
  * @param[out] S CIO locator [rad]. */
  virtual void earthOrientationParameter(const Time &timeGPS, Double &xp, Double &yp, Double &sp, Double &deltaUT, Double &LOD, Double &X, Double &Y, Double &S) const;

  /** @brief creates an derived instance of this class. */
  static EarthRotationPtr create(Config &config, const std::string &name);
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class EarthRotation.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a earthRotation is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] earthRotation Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates EarthRotation */
template<> Bool readConfig(Config &config, const std::string &name, EarthRotationPtr &earthRotation, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***********************************************/

#endif
