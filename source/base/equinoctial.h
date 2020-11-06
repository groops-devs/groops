/***********************************************/
/**
* @file equinoctial.h
*
* @brief Equinoctial elements and orbits.
*
* @author Matthias Ellmer
* @date 2015-06-10
*
* @see Broucke, R. A., and P. J. Cefola. 1972. “On the Equinoctial Orbit Elements.” Celestial Mechanics 5 (3): 303–10. doi:10.1007/BF01228432.
* @see Broucke, R. A. 1970. “On the Matrizant of the Two-Body Problem.” Astronomy and Astrophysics 6 (June): 173–82.
* @see Montenbruck, Oliver, and Eberhard Gill. 2000. Satellite Orbits. Springer-Verlag Berlin Heidelberg New York. Chapter 2.2.5, pp. 30-32
* @see Danielson, Donald A., Christopher Patrick Sagovac, Beny Neta, and Leo W. Early. 1995. “Semianalytic Satellite Theory.” DTIC Document.
*/
/***********************************************/

#ifndef __GROOPS_EQUINOCTIAL__
#define __GROOPS_EQUINOCTIAL__

#include "base/importStd.h"
#include "base/time.h"
#include "base/vector3d.h"

class Kepler;

/***** CLASS ***********************************/

/** @brief Equinoctial orbits.
* @ingroup base */
class Equinoctial
{
public:
  Double GM;      //!< Specific gravitational constant
  Double a;       //!< a                     | Semi-major axis
  Double h;       //!< e sin(omega + Omega)
  Double k;       //!< e cos(omega + Omega)
  Double p;       //!< tan(i/2) + sin(Omega)
  Double q;       //!< tan(i/2) + cos(Omega)
  LongDouble l0;  //!< M0 + omega + Omega    | Mean longitude at t = tau
  Time   tau;     //!<                       | Time of perigee

  /// Default constructor.
  Equinoctial();

  /** @brief Constructor with equinoctial elements
  * @param a    Semi-major axis
  * @param h
  * @param k
  * @param p
  * @param q
  * @param l0   Mean longitude at time of perigee
  * @param tau  Time of perigee
  * @param GM   Specific gravitational constant */
 Equinoctial(Double a, Double h, Double k, Double p, Double q, Double l0, const Time &tau, Double GM = DEFAULT_GM);

  /** @brief Equinoctial from Kepler orbit
  * @param kepler Kepler object */
  explicit Equinoctial(const Kepler &kepler);

  /** @brief Constructor with state vector at initial time @a epoch.
  * @param epoch    Reference epoch. Position and velocity are given at this epoch. The motion will be parametrized with respect to this epoch.
  * @param position Initial position
  * @param velocity Initial velocity
  * @param GM   Specific gravitational constant */
  explicit Equinoctial(const Time &epoch, const Vector3d &position, const Vector3d &velocity, Double GM=DEFAULT_GM);

  /** @brief Constructor with state vector at arbitrary time @a time.
  * @param epoch    Reference epoch. The motion will be parametrized with respect to this epoch.
  * @param time     Time of position and velocity
  * @param position Initial position
  * @param velocity Initial velocity
  * @param GM   Specific gravitational constant
  * @see Chapter 2.1.5 Danielson, Donald A., Christopher Patrick Sagovac, Beny Neta, and Leo W. Early. 1995. “Semianalytic Satellite Theory.” */
  explicit Equinoctial(const Time &epoch, const Time &time, const Vector3d &position, const Vector3d &velocity, Double GM=DEFAULT_GM);

  /** @brief Mean longitude @a lambda at @a time
  * @param time Time of mean longitude
  * @see Eq. 14.a, Broucke, R. A., and P. J. Cefola. 1972. “On the Equinoctial Orbit Elements.” Celestial Mechanics 5 (3): 303–10. doi:10.1007/BF01228432.
  * We use a LongDouble for numerical reasons. Roundoff in time differences can lead to discrepancies in
  * computation of different orbits. We use LongDouble for the mean Longitude to get around this. Needs to be considered here
  * and in constructor from position and velocity. */
  LongDouble meanLongitude(const Time &time) const;

  /** @brief Eccentric longitude @a F at @a time through Newton's method
  * @param time Eccentric longitude at this time
  * @see Eq. 14.b, Broucke, R. A., and P. J. Cefola. 1972. “On the Equinoctial Orbit Elements.” Celestial Mechanics 5 (3): 303–10. doi:10.1007/BF01228432.
  */
  Double eccentricLongitude(const Time &time) const;

  /** @brief True longitude @a L at @a time
  * @param time True longitude at this time
  * @see Eq. 14.c, Broucke, R. A., and P. J. Cefola. 1972. “On the Equinoctial Orbit Elements.” Celestial Mechanics 5 (3): 303–10. doi:10.1007/BF01228432.
  * @see Chapter 2.1.4 Eq 5, Danielson, Donald A., Christopher Patrick Sagovac, Beny Neta, and Leo W. Early. 1995. “Semianalytic Satellite Theory.” DTIC Document */
  Double trueLongitude(const Time &time) const;

  /** @brief Mean anomaly @a M at @a time
  * @param time Time of mean anomaly
  * @see Chapter 2.1.3 Eq 2, Danielson, Donald A., Christopher Patrick Sagovac, Beny Neta, and Leo W. Early. 1995. “Semianalytic Satellite Theory.” DTIC Document */
  Double meanAnomaly(const Time &time) const;

  /** @brief Orthonormal vectors defining the orbital plane system
  * @param[out] F perigee direction
  * @param[out] G corresponding to true anomaly 90°
  * @param[out] H completing the system
  * @see Eq. 2.71 ff in Montenbruck, Oliver, and Eberhard Gill. 2000. Satellite Orbits. Springer-Verlag Berlin Heidelberg New York. */
  void orthonormals(Vector3d &F, Vector3d &G, Vector3d &H) const;

  /** @brief In-plane orthonormal vectors
  * @param[out] F perigee direction
  * @param[out] G corresponding to true anomaly 90°
  * @see Eq. 2.71 ff in Montenbruck, Oliver, and Eberhard Gill. 2000. Satellite Orbits. Springer-Verlag Berlin Heidelberg New York. */
  void orthonormals(Vector3d &F, Vector3d &G) const;

  /** @brief State transition matrix and inverse at @a time.
  * @param      time   Evaluation time
  * @returns           State transition matrix. The rows contain the partial derivatives of position and velocity with regard to equinoctial elements.
  * @see Chapter 2.1.6 and 2.1.7 -- Danielson, Donald A., Christopher Patrick Sagovac, Beny Neta, and Leo W. Early. 1995. “Semianalytic Satellite Theory.” DTIC Document. */
  Matrix stateTransitionMatrix(const Time &time);

  /** @brief Position at @a time.
  * @param time time
  * @returns position at @a time */
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
  * @param[out] acceleration acceleration at @a time
  * @see Eq. 2.71 in Montenbruck, Oliver, and Eberhard Gill. 2000. Satellite Orbits. Springer-Verlag Berlin Heidelberg New York. */
  void orbit(const Time &time, Vector3d &position, Vector3d &velocity, Vector3d &acceleration) const;

  /** @brief Position and velocity in orbital plane at @a time
  * @param time time
  * @param[out] position position at @a time
  * @param[out] velocity velocity at @a time
  * @see Eq. 2.71 in Montenbruck, Oliver, and Eberhard Gill. 2000. Satellite Orbits. Springer-Verlag Berlin Heidelberg New York. */
  void inPlaneOrbit(const Time &time, Vector3d &position, Vector3d &velocity) const;

  /** @brief Position, velocity, and acceleration in orbital plane at @a time
  * @param time time
  * @param[out] position position at @a time
  * @param[out] velocity velocity at @a time
  * @param[out] acceleration acceleration at @a time
  * @see Eq. 2.71 in Montenbruck, Oliver, and Eberhard Gill. 2000. Satellite Orbits. Springer-Verlag Berlin Heidelberg New York. */
  void inPlaneOrbit(const Time &time, Vector3d &position, Vector3d &velocity, Vector3d &acceleration) const;
};

/***********************************************/

#endif /* __GROOPS_KEPLER__ */
