/***********************************************/
/**
* @file forces.h
*
* @brief Wraps forces and force-generating potentials.
*
* @author Matthias Ellmer
* @date 2017-01-18
*
*/
/***********************************************/

#ifndef __GROOPS_FORCES__
#define __GROOPS_FORCES__

// Latex documentation
#ifdef DOCSTRING_Forces
static const char *docstringForces = R"(
\section{Forces}\label{forcesType}
This class provides the forces acting on a satellite.
This encompasses \configClass{gravityfield}{gravityfieldType}, \configClass{tides}{tidesType}
and \configClass{miscAccelerations}{miscAccelerationsType}.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "files/fileSatelliteModel.h"
#include "classes/ephemerides/ephemerides.h"
#include "classes/gravityfield/gravityfield.h"
#include "classes/miscAccelerations/miscAccelerations.h"
#include "classes/tides/tides.h"

/**
* @defgroup forcesGroup Forces
* @brief Forces acting on satellites.
* @ingroup classesGroup
* The interface is given by @ref Forces.
* An Instance can be created by @ref readConfig. */
/// @{

/***** TYPES ***********************************/

class Forces;
typedef std::shared_ptr<Forces> ForcesPtr;

/***** CLASS ***********************************/

/** @brief Wrapper around forces and force-generating types.
* An Instance of this class can be created by @ref readConfig. */
class Forces
{
public:
  /// Constructor.
  Forces(Config &config, const std::string &name);

  /** @brief Compute full acceleration in TRF
  * @param satellite model for misc accelerations
  * @param time Time.
  * @param position in CRF [m].
  * @param velocity in CRF [m/s].
  * @param rotSat   Sat -> CRF
  * @param rotEarth CRF -> TRF
  * @param rotation need for computation of polar motion.
  * @param ephemerides Position of Sun and Moon.
  * @return acceleration in TRF(!) [m/s^2] */
  Vector3d acceleration(SatelliteModelPtr satellite, const Time &time, const Vector3d &position, const Vector3d &velocity,
                        const Rotary3d &rotSat, const Rotary3d &rotEarth, EarthRotationPtr rotation, EphemeridesPtr ephemerides) const;

  /** @brief creates an derived instance of this class. */
  static ForcesPtr create(Config &config, const std::string &name) {return ForcesPtr(new Forces(config, name));}

private:
  GravityfieldPtr      gravityfield;
  TidesPtr             tides;
  MiscAccelerationsPtr miscAccelerations;
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class Forces.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a forces is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] forces Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates Forces */
template<> Bool readConfig(Config &config, const std::string &name, ForcesPtr &forces, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***********************************************/

#endif // __GROOPS_FORCES__

