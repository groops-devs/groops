/***********************************************/
/**
* @file miscAccelerations.h
*
* @brief Non conservative forces acting on satellites.
*
* @author Torsten Mayer-Guerr
* @date 2013-02-09
*
*/
/***********************************************/

#ifndef __GROOPS_MISCACCELERATIONS__
#define __GROOPS_MISCACCELERATIONS__

// Latex documentation
#ifdef DOCSTRING_MiscAccelerations
static const char *docstringMiscAccelerations = R"(
\section{MiscAccelerations}\label{miscAccelerationsType}
This class gives the non conservative forces acting on satellites.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "files/fileSatelliteModel.h"
#include "classes/ephemerides/ephemerides.h"

/**
* @defgroup miscAccelerationsGroup MiscAccelerations
* @brief Non conservative forces acting on satellites.
* @ingroup classesGroup
* The interface is given by @ref MiscAccelerations.
* An Instance can be created by @ref readConfig. */
/// @{

/***** TYPES ***********************************/

class MiscAccelerations;
class MiscAccelerationsBase;
typedef std::shared_ptr<MiscAccelerations> MiscAccelerationsPtr;

/***** CLASS ***********************************/

/** @brief Non conservative forces acting on satellites.
* An Instance of this class can be created by @ref readConfig. */
class MiscAccelerations
{
  std::vector<MiscAccelerationsBase*> acc;

public:
  /// Constructor.
  MiscAccelerations(Config &config, const std::string &name);

  /// Destructor.
  ~MiscAccelerations();

  /** @brief Acceleration.
  * @param satellite Macro model.
  * @param time Time.
  * @param position in CRF [m].
  * @param velocity in CRF [m/s].
  * @param rotSat   Sat -> CRF
  * @param rotEarth CRF -> TRF
  * @param ephemerides Position of Sun and Moon.
  * @return acceleration in TRF(!) [m/s^2] */
  Vector3d acceleration(SatelliteModelPtr satellite, const Time &time, const Vector3d &position, const Vector3d &velocity,
                        const Rotary3d &rotSat, const Rotary3d &rotEarth, EphemeridesPtr ephemerides);

  /** @brief creates an derived instance of this class. */
  static MiscAccelerationsPtr create(Config &config, const std::string &name) {return MiscAccelerationsPtr(new MiscAccelerations(config, name));}
};


/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class MiscAccelerations.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a acceleration without forces is created.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] acceleration Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates MiscAccelerations */
template<> Bool readConfig(Config &config, const std::string &name, MiscAccelerationsPtr &acceleration, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***** CLASS ***********************************/

// Internal class
class MiscAccelerationsBase
{
public:
virtual ~MiscAccelerationsBase() {}
virtual  Vector3d acceleration(SatelliteModelPtr satellite, const Time &time, const Vector3d &position, const Vector3d &velocity,
                               const Rotary3d &rotSat, const Rotary3d &rotEarth, EphemeridesPtr ephemerides) = 0;
};

/***********************************************/

#endif
