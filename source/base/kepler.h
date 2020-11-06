/***********************************************/
/**
* @file kepler.h
*
* @brief Keplerian elements and orbits.
*
* @author Torsten Mayer-Guerr
* @date 2004-09-12
*
*/
/***********************************************/

#ifndef __GROOPS_KEPLER__
#define __GROOPS_KEPLER__

#include "base/importStd.h"
#include "base/vector3d.h"
#include "base/time.h"

class Equinoctial;

/***** CLASS ***********************************/

/** @brief Keplerian orbits.
* @ingroup base */
class Kepler
{
public:
  Double GM;    //!< Specific gravitational constant
  Double Omega; //!< ascending node
  Double i;     //!< inclination
  Double omega; //!< argument of perigee
  Double a;     //!< major semi axis
  Double e;     //!< excentricity
  Time   tau;   //!< time of perigee

  /// Default constructor.
  Kepler();

  /** @brief Constructor with keplerian elements.
  * @param Omega ascending node
  * @param i     inclination
  * @param omega argument of perigee
  * @param a     major semi axis
  * @param e     excentricity
  * @param tau   time of perigee
  * @param GM    specific gravitational constant */
  Kepler(Double Omega, Double i, Double omega, Double a, Double e, const Time &tau, Double GM=DEFAULT_GM);

  /** @brief Constructor with keplerian elements.
  * @param time  time of mean anomaly
  * @param Omega ascending node
  * @param i     inclination
  * @param omega argument of perigee
  * @param a     major semi axis
  * @param e     excentricity
  * @param M     mean anomaly
  * @param GM    specific gravitational constant */
  Kepler(const Time &time, Double Omega, Double i, Double omega, Double a, Double e, Double M, Double GM=DEFAULT_GM);

  /** @brief Constructor with initial state vector.
  * @param time Start time
  * @param position Start position
  * @param velocity Start velocity
  * @param GM   Specific gravitational constant */
  Kepler(const Time &time, const Vector3d &position, const Vector3d &velocity, Double GM=DEFAULT_GM);

  /** @brief Kepler from Equinoctial orbit
  * @param equinoctial Equinoctial object
  * @see Chapter 2.1.3 Danielson, Donald A., Christopher Patrick Sagovac, Beny Neta, and Leo W. Early. 1995. “Semianalytic Satellite Theory.” */
  explicit Kepler(const Equinoctial &equinoctial);

  /// Mean anomaly at @a time.
  Double meanAnomaly(const Time &time) const;

  /// True anomaly at @a time.
  Double trueAnomaly(const Time &time) const;

  /// Position at @a time.
  const Vector3d position(const Time &time) const;

  /** @brief Position and velocity at @a time.
  * @param time time
  * @param[out] position position at @a time
  * @param[out] velocity velocity at @a time */
  void orbit(const Time &time, Vector3d &position, Vector3d &velocity) const;

  /** @brief Position, velocity, and acceleration at @a time.
  * @param time time
  * @param[out] position position at @a time
  * @param[out] velocity velocity at @a time
  * @param[out] acceleration acceleration at @a time */
  void orbit(const Time &time, Vector3d &position, Vector3d &velocity, Vector3d &acceleration) const;
};

/***********************************************/

#endif /* __GROOPS_KEPLER__ */
