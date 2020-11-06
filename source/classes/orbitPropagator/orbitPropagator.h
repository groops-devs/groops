/***********************************************/
/**
* @file orbitPropagator.h
*
* @brief Propagate orbit using different algorithms.
*
* @author Matthias Ellmer
* @date 2017-01-19
*
*/
/***********************************************/

#ifndef __GROOPS_ORBITPROPAGATOR__
#define __GROOPS_ORBITPROPAGATOR__

// Latex documentation
#ifdef DOCSTRING_OrbitPropagator
static const char *docstringOrbitPropagator = R"(
\section{OrbitPropagator}\label{orbitPropagatorType}
Implements the propagation of a satellite orbit under
the influence of \configClass{forces}{forcesType} as
used in \program{SimulateOrbit}
(dynamic orbits from numerical orbit integration).
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "files/fileInstrument.h"
#include "files/fileSatelliteModel.h"
#include "classes/forces/forces.h"

/**
* @defgroup orbitPropagatorGroup OrbitPropagator
* @brief Propagate a satellite orbit.
* @ingroup classesGroup
* The interface is given by @ref OrbitPropagator.
* An Instance can be created by @ref readConfig. */
/// @{

/***** TYPES ***********************************/

class OrbitPropagator;
typedef std::shared_ptr<OrbitPropagator> OrbitPropagatorPtr;

/***** CLASS ***********************************/

/** @brief Propagate a satellite orbit.
* An Instance of this class can be created by @ref readConfig. */
class OrbitPropagator
{
public:
  /// Destructor.
  virtual ~OrbitPropagator() {}

  /** @brief Propagate orbit arc through force field
  * @param startEpoch Initial epoch of the orbit
  * @param sampling Time difference between epochs
  * @param posCount Total number of epochs in arc (including @a startEpoch)
  * @param forces Force models to consider
  * @param satellite Satellite macro model to consider
  * @param earthRotation for rotations in forces evaluation
  * @param ephemerides Position of Sun and Moon.
  * @param timing Show a timer.
  * @returns OrbitArc arc of size @a posCount, with @a startEpoch at index 0  */
  virtual OrbitArc integrateArc(OrbitEpoch startEpoch, Time sampling, UInt posCount, ForcesPtr forces, SatelliteModelPtr satellite,
                                EarthRotationPtr earthRotation, EphemeridesPtr ephemerides, Bool timing=TRUE) const = 0;

  /** @brief Generic orientation definition in orbital system.
  * @param time of evaluation
  * @param position of the satellite
  * @param velocity of the satellite
  * @param satellite model of the satellite
  * @returns Rotary3d satellite orientation, with x as along-track, y cross-track and z perpendicular.  */
  virtual Rotary3d orientation(const Time &time, const Vector3d &position, const Vector3d &velocity, SatelliteModelPtr satellite) const;

  static Arc flip(const Arc &arc);

  /** @brief creates a derived instance of this class. */
  static OrbitPropagatorPtr create(Config &config, const std::string &name);

protected:
  Vector3d acceleration(const OrbitEpoch &epoch, ForcesPtr forces, SatelliteModelPtr satellite, EarthRotationPtr earthRotation, EphemeridesPtr ephemerides) const;
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class OrbitPropagator.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a orbitPropagator is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] orbitPropagator Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates OrbitPropagator */
template<> Bool readConfig(Config &config, const std::string &name, OrbitPropagatorPtr &orbitPropagator, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

#endif
