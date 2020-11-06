/***********************************************/
/**
* @file constants.h
*
* @brief define constants.
*
* Constants are only valid after call of main().
* New constants must be set in inputOutput/settings.cpp.
*
* @author Torsten Mayer-Guerr
* @date 2010-03-01
*
*/
/***********************************************/

#ifndef __GROOPS_CONSTANTS__
#define __GROOPS_CONSTANTS__

#include "base/importStd.h"

/**
* @defgroup constants Constants
* @brief General constants used in groops.
* @ingroup base
* Constants can be overwritten with groopsDefaults.xml file on input. */
// @{

/***** CONST ***********************************/

constexpr Double PI             = 3.141592653589793238462643383279502884; //!< 3.1415.
constexpr Double RAD2DEG        = 180.0/PI;           //!< Conversion radian -> degree
constexpr Double DEG2RAD        = PI/180.0;           //!< Conversion degree -> radian.
constexpr UInt   MAX_UINT       = 4294967295U;        //!< max. representable integer number.
constexpr UInt   INFINITYDEGREE = MAX_UINT;           //!< infinity polynomial/harmonics degree.
constexpr UInt   NULLINDEX      = MAX_UINT;           //!< undefined index.
constexpr Double NAN_EXPR       = std::numeric_limits<Double>::has_quiet_NaN ? std::numeric_limits<Double>::quiet_NaN() : NAN; //!< not a number symbol

[[deprecated("Use RAD2DEG instead (or 1/RHO -> RAD2DEG")]] constexpr Double RHO = RAD2DEG;

extern Double LIGHT_VELOCITY;  //!< speed of light [m/s].

/***** DEFAULTS ********************************/

extern std::string STRING_DEFAULT_GM;       //!< Gravitational constant of the Earth (EGM96).
extern std::string STRING_DEFAULT_R;        //!< Reference radius (EGM96).
extern std::string STRING_DEFAULT_GRS80_a;  //!< Semi major axis of Earth's ellipsoid.
extern std::string STRING_DEFAULT_GRS80_f;  //!< inverse flattening of Earth's ellipsoid.

extern Double DEFAULT_GM;      //!< Gravitational constant of the Earth (EGM96).
extern Double DEFAULT_R;       //!< Reference radius (EGM96).
extern Double DEFAULT_GRS80_a; //!< Semi major axis of Earth's ellipsoid.
extern Double DEFAULT_GRS80_f; //!< inverse flattening of Earth's ellipsoid.

/***** CONSTANTS for Planets and Tides *********/

extern Double GRAVITATIONALCONSTANT; //!< Gravitational constant 6.673e-11.
extern Double R_Earth;               //!< Earth radius (Used in Planets).
extern Double R_Moon;                //!< Moon radius.
extern Double GM_Earth;              //!< Gravitational constant of the Earth (Used in Planets).
extern Double GM_Sun;                //!< Gravitational constant of the sun (Used in Planets).
extern Double GM_Moon;               //!< Gravitational constant of the moon (Used in Planets).
extern Double GM_MERCURY;            //!< Gravitational constant of mercury (Used in Planets).
extern Double GM_VENUS;              //!< Gravitational constant of venus (Used in Planets).
extern Double GM_MARS;               //!< Gravitational constant of mars (Used in Planets).
extern Double GM_JUPITER;            //!< Gravitational constant of jupiter (Used in Planets).
extern Double GM_SATURN;             //!< Gravitational constant of saturn (Used in Planets).

/***** Ionosphere ******************************/

namespace Ionosphere
{
  constexpr Double e0 =  8.854187817e-12;              //!< Permittivity of free space [F/m].
  constexpr Double me =  9.1093835611e-31;             //!< Electron mass [kg]
  constexpr Double E  = -1.602176634e-19;              //!< Electron charge [C]
  constexpr Double Ap =  1e16*(E*E/(4*PI*PI*e0*me))/2; //!< 40.3e16 in [TECU -> m]
}

/***** TIME SHIFTS *****************************/

extern Double      TIME_EPSILON;   //!< Margin to consider two times identical [seconds].
extern Double      DELTA_TAI_GPS;  //!< Seconds to be added to GPS time to get TAI.
extern Double      DELTA_TT_GPS;   //!< Seconds to be added to GPS time to get TT.
extern Double      J2000;          //!< J2000 in mjd
extern std::string STRING_J2000;   //!< J2000 in mjd

extern std::vector<Int>    MJD_UTC_GPS;   //!< Time of new introduced leap second.
extern std::vector<Double> DELTA_UTC_GPS; //!< Leap seconds.

// @} group constants

/***********************************************/

#endif /* __GROOPS__ */
