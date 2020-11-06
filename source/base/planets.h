/***********************************************/
/**
* @file planets.h
*
* @brief Planet and Earth fundamentals.
*
* @author Torsten Mayer-Guerr
* @date 2005-04-26
*
*/
/***********************************************/

#ifndef __GROOPS_PLANETS__
#define __GROOPS_PLANETS__

#include "base/importStd.h"
#include "base/matrix.h"
#include "base/time.h"
#include "base/vector3d.h"
#include "base/rotary3d.h"

class FileName;

/***** CLASS ***********************************/

/** @brief planet and Earth fundamentals.
* @ingroup base */
namespace Planets
{
  /** @brief Fundamental arguments of nutation theory.
  * -  F(0) = l  = Mean Anomaly of the Moon
  * -  F(1) = l' = Mean Anomaly of the Sun
  * -  F(2) = L-Omega
  * -  F(3) = D = Mean Elongation of the Moon from Sun
  * -  F(4) = Omega = Mean Longitude of the Ascening Node of the Moon */
  Vector fundamentals(const Time &timeGPS);

  /** @brief Position of sun in celestial frame.
  * Reference:  The Astronomical Almanac, page C24.
  * Accurate to 0.01 degrees between the years 1950 through 2050. */
  Vector3d positionSun(const Time &timeGPS);

  /** @brief Greenwich Mean Sideral Time.
  * roation angle (CRF -> TRF) [rad]. */
  Double gmst(const Time &timeUT1);

  /** @brief Earth rotation angle.
  * roation angle (CRF -> TRF) [rad]. */
  Double ERA(const Time &timeUT1);

  /** @brief Approximate Earth rotation.
  * Without Earth rotation parameter (EOP). */
  Rotary3d celestial2TerrestrialFrame(const Time &timeGPS);

  /** @brief Gravity of the Earth.
  * Normal Gravity Formula (Geodetic Reference System 1980). */
  Double normalGravity(const Vector3d &p);
} // end Planets

/***********************************************/

#endif /* __GROOPS_PLANET__ */
