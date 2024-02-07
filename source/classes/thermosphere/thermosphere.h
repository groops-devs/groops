/***********************************************/
/**
* @file thermosphere.h
*
* @brief Density, temperature and velocity.
*
* @author Torsten Mayer-Guerr
* @date 2020-02-20
*
*/
/***********************************************/

#ifndef __GROOPS_THERMOSPHERE__
#define __GROOPS_THERMOSPHERE__

// Latex documentation
#ifdef DOCSTRING_Thermosphere
static const char *docstringThermosphere = R"(
\section{Thermosphere}\label{thermosphereType}
This class provides functions for calculating the density, temperature and velocity
in the thermosphere.
The wind is computed by HWM14 model if \config{hwm14DataDirectory} is provided.
A quiet thermosphere is assumed if \config{inputfileMagnetic3hAp} is not given.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "files/fileInstrument.h"

/**
* @defgroup thermosphereGroup Thermosphere
* @brief Density, temperature and velocity.
* @ingroup classesGroup
* The interface is given by @ref Thermosphere. */
/// @{

/***** TYPES ***********************************/

class Thermosphere;
typedef std::shared_ptr<Thermosphere> ThermospherePtr;

/***** CLASS ***********************************/

/** @brief Density, temperature and velocity.
* An Instance of this class can be created by @ref readConfig. */
class Thermosphere
{
public:
  /// Destructor.
  virtual ~Thermosphere() {}

  /** @brief Thermospheric state.
  * @param time GPS time
  * @param position in TRF [m]
  * @param[out] density  [kg/m^3]
  * @param[out] temperature  [K]
  * @param[out] velocity wind in TRF [m/s] */
  virtual void state(const Time &time, const Vector3d &position, Double &density, Double &temperature, Vector3d &velocity) const = 0;

  /** @brief creates an derived instance of this class. */
  static ThermospherePtr create(Config &config, const std::string &name);

protected:
  FileName fileNameHwm14Path;

  static Vector getIndices(const MiscValuesArc &arc, const Time &time, Bool interpolate);

  MiscValuesArc magnetic3hAp;
  Vector3d wind(const Time &time, const Vector3d &position) const;
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class Thermosphere.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a thermosphere is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] thermosphere Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates Thermosphere */
template<> Bool readConfig(Config &config, const std::string &name, ThermospherePtr &thermosphere, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***********************************************/

#endif /* __GROOPS__ */
