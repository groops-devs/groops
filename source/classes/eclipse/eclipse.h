/***********************************************/
/**
* @file eclipse.h
*
* @brief Shadowing of satellites by moon and Earth.
*
* @author Torsten Mayer-Guerr
* @date 2020-03-08
*
*/
/***********************************************/

#ifndef __GROOPS_ECLIPSE__
#define __GROOPS_ECLIPSE__

// Latex documentation
#ifdef DOCSTRING_Eclipse
static const char *docstringEclipse = R"(
\section{Eclipse}\label{eclipseType}
Shadowing of satellites by moon and Earth provided as factor
between $[0,1]$ with 0: full shadow and 1: full sun light.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "classes/ephemerides/ephemerides.h"

/**
* @defgroup eclipseGroup Eclipse
* @brief Shadowing of satellites by moon and Earth.
* @ingroup classesGroup
* The interface is given by @ref Eclipse.
* An Instance can be created by @ref readConfig. */
/// @{

/***** TYPES ***********************************/

class Eclipse;
typedef std::shared_ptr<Eclipse> EclipsePtr;

/***** CLASS ***********************************/

/** @brief Shadowing of satellites by moon and Earth.
* An Instance of this class can be created by @ref readConfig. */
class Eclipse
{
public:
  /// Destructor.
  virtual ~Eclipse() {}

  /** @brief Scaling factor for satellite during (Earth/Moon) shadow crossing.
  * @param timeGPS  Modified julian date (MJD) in GPS time system.
  * @param position Position of the satellite in celestial frame (CRF).
  * @param ephemerides Position of Sun and Moon.
  * @return Scaling factor [0..1] (0 = full shadow, 1 = no shadow, between 0 and 1 = partial shadow) */
  virtual Double factor(const Time &timeGPS, const Vector3d &position, EphemeridesPtr ephemerides) const = 0;

  /** @brief creates an derived instance of this class. */
  static EclipsePtr create(Config &config, const std::string &name);

protected:
  static Double shadowScalingFactor(const Vector3d &posSat, const Vector3d &posSun, const Vector3d &posBody, const Double &radiusBody);
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class Eclipse.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a eclipse is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] eclipse Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates Eclipse */
template<> Bool readConfig(Config &config, const std::string &name, EclipsePtr &eclipse, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***********************************************/

#endif
