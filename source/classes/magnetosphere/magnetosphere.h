/***********************************************/
/**
* @file magnetosphere.h
*
* @brief Magentic field of the Earth.
*
* @author Torsten Mayer-Guerr
* @date 2020-06-08
*
*/
/***********************************************/

#ifndef __GROOPS_MAGNETOSPHERE__
#define __GROOPS_MAGNETOSPHERE__

// Latex documentation
#ifdef DOCSTRING_Magnetosphere
static const char *docstringMagnetosphere = R"(
\section{Magnetosphere}\label{magnetosphereType}
This class provides functions of the magnetic field of the Earth.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"

/**
* @defgroup magnetosphereGroup Magnetosphere
* @brief Magentic field of the Earth.
* @ingroup classesGroup
* The interface is given by @ref Magnetosphere.
* An Instance can be created by @ref readConfig. */
/// @{

/***** TYPES ***********************************/

class Magnetosphere;
typedef std::shared_ptr<Magnetosphere> MagnetospherePtr;

/***** CLASS ***********************************/

/** @brief Magentic field of the Earth.
* An Instance of this class can be created by @ref readConfig. */
class Magnetosphere
{
public:
  /// Destructor.
  virtual ~Magnetosphere() {}

  /** @brief Magnetic field of the Earth.
  * @param time GPS time
  * @param position in TRF [m]
  * @return Magentic vector in terrestrial reference system (TRF) [Tesla = kg/A/s^2]. */
  virtual Vector3d magenticFieldVector(const Time &time, const Vector3d &position) const = 0;

  /** @brief Geomagnetic north pole in terrestrial frame (TRF).
  * Unit vector. */
  virtual Vector3d geomagneticNorthPole(const Time &time) const = 0;

  /** @brief Rotation from celestial frame (CRF) to solar geomagnetic frame (SGF). */
  virtual Rotary3d rotaryCelestial2SolarGeomagneticFrame(const Time &time) const;

  /** @brief creates an derived instance of this class. */
  static MagnetospherePtr create(Config &config, const std::string &name);
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class Magnetosphere.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a magnetosphere is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] magnetosphere Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates Magnetosphere */
template<> Bool readConfig(Config &config, const std::string &name, MagnetospherePtr &magnetosphere, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***********************************************/

#endif /* __GROOPS__ */
