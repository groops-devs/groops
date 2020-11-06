/***********************************************/
/**
* @file ephemerides.h
*
* @brief Ephemerides of Sun, Moon and planets.
*
* @author Torsten Mayer-Guerr
* @date 2020-06-07
*
*/
/***********************************************/

#ifndef __GROOPS_EPHEMERIDES__
#define __GROOPS_EPHEMERIDES__

// Latex documentation
#ifdef DOCSTRING_Ephemerides
static const char *docstringEphemerides = R"(
\section{Ephemerides}\label{ephemeridesType}
Ephemerides of Sun, Moon and planets.
The coordinate system is defined as center of \configClass{origin}{planetType}.
)";
#endif

// Latex documentation
#ifdef DOCSTRING_EphemeridesPlanet
static const char *docstringEphemeridesPlanet = R"(
\section{Planet}\label{planetType}
Defines the planet to compute the \configClass{ephemeris}{ephemeridesType}.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"

/**
* @defgroup ephemeridesGroup Ephemerides
* @brief Ephemerides of Sun, Moon and planets.
* @ingroup classesGroup
* The interface is given by @ref Ephemerides.
* An Instance can be created by @ref readConfig. */
/// @{

/***** TYPES ***********************************/

class Ephemerides;
typedef std::shared_ptr<Ephemerides> EphemeridesPtr;

/***** CLASS ***********************************/

/** @brief Ephemerides of Sun, Moon and planets.
* An Instance of this class can be created by @ref readConfig. */
class Ephemerides
{
public:
  /// planet identifier
  enum Planet {MERCURY             = 1,
               VENUS               = 2,
               EARTH               = 3,
               MARS                = 4,
               JUPITER             = 5,
               SATURN              = 6,
               URANUS              = 7,
               NEPTUNE             = 8,
               PLUTO               = 9,
               MOON                = 10,
               SUN                 = 11,
               SOLARBARYCENTER     = 12,
               EARTHMOONBARYCENTER = 13};

  /// Destructor.
  virtual ~Ephemerides() {}

  /** @brief Position of a planet in celestial reference system (CRF).  */
  virtual Vector3d position(const Time &timeGPS, Planet planet) = 0;

  /** @brief Position and velocity of a planet in celestial reference system (CRF).  */
  virtual void ephemeris(const Time &timeGPS, Planet planet, Vector3d &position, Vector3d &velocity) = 0;

  virtual Planet origin() const = 0;

  /** @brief creates an derived instance of this class. */
  static EphemeridesPtr create(Config &config, const std::string &name);
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class Ephemerides.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a var is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] var Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates Ephemerides */
template<> Bool readConfig(Config &config, const std::string &name, EphemeridesPtr &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/** @brief Reads an planet from config.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a var is untouched.
* @param config The config node which includes the node with the options for @a var
* @param name Tag name in the config.
* @param[out] var read from config.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Used for @a var if no node is in config.
* @param annotation Short description. */
template<> Bool readConfig(Config &config, const std::string &name, Ephemerides::Planet &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***********************************************/

#endif /* __GROOPS__ */
