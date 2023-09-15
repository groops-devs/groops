/***********************************************/
/**
* @file orbitPropagatorRungeKutta4.h
*
* @brief Propagate a dynamic orbit using the classical Runge-Kutta 4 method.
* @see orbitPropagator
*
* @author Matthias Ellmer
* @date 2017-01-19
*
*/
/***********************************************/

#ifndef __GROOPS_ORBITPROPAGATORRUNGEKUTTA4__
#define __GROOPS_ORBITPROPAGATORRUNGEKUTTA4__

// Latex documentation
#ifdef DOCSTRING_OrbitPropagator
static const char *docstringOrbitPropagatorRungeKutta4 = R"(
\subsection{RungeKutta4}
This class implements the classical Runge-Kutta 4 method of orbit propagation
for satellite orbit under the influence of \configClass{Forces}{forcesType}.
No step-width control or other advanced features are implemented.
Satellite is assumed to be oriented along-track.
See: Montenbruck, Oliver, and Eberhard Gill. 2000. Satellite Orbits
)";
#endif

/***********************************************/

#include "classes/orbitPropagator/orbitPropagator.h"

/***** CLASS ***********************************/

/** @brief Propagate orbit using the Runge-Kutta 4 method.
* @ingroup orbitPropagatorGroup
* @see orbitPropagator */
class OrbitPropagatorRungeKutta4 : public OrbitPropagator
{
public:
  OrbitPropagatorRungeKutta4(Config &/*config*/) {}
  OrbitArc integrateArc(const OrbitEpoch &startEpoch, const Time &sampling, UInt posCount, ForcesPtr forces, SatelliteModelPtr satellite,
                        EarthRotationPtr earthRotation, EphemeridesPtr ephemerides, Bool timing) const override;
};

/***********************************************/

inline OrbitArc OrbitPropagatorRungeKutta4::integrateArc(const OrbitEpoch &startEpoch, const Time &sampling, UInt posCount, ForcesPtr forces,
                                                         SatelliteModelPtr satellite, EarthRotationPtr earthRotation, EphemeridesPtr ephemerides, Bool timing) const
{
  try
  {
    OrbitArc orbit;
    orbit.push_back(startEpoch);
    orbit.back().acceleration = acceleration(startEpoch, forces, satellite, earthRotation, ephemerides);
    const Double dt = sampling.seconds();

    Single::forEach(posCount-1, [&](UInt k)
    {
      // Evaluate accelerations at 4 positions between current epoch and next
      OrbitEpoch k1 = orbit.back();

      OrbitEpoch k2 = k1;
      k2.time        += seconds2time(dt/2.);
      k2.position    += dt/2. * k1.velocity;
      k2.velocity    += dt/2. * k1.acceleration;
      k2.acceleration = acceleration(k2, forces, satellite, earthRotation, ephemerides);

      OrbitEpoch k3 = k1;
      k3.time        += seconds2time(dt/2.);
      k3.position    += dt/2. * k2.velocity;
      k3.velocity    += dt/2. * k2.acceleration;
      k3.acceleration = acceleration(k3, forces, satellite, earthRotation, ephemerides);

      OrbitEpoch k4 = k1;
      k4.time        += sampling;
      k4.position    += dt * k3.velocity;
      k4.velocity    += dt * k3.acceleration;
      k4.acceleration = acceleration(k4, forces, satellite, earthRotation, ephemerides);

      // Compute final value for this epoch
      OrbitEpoch epoch = k1;
      epoch.time         = startEpoch.time + (k+1)*sampling;
      epoch.position    += (dt/6.) * (k1.velocity + 2*k2.velocity + 2*k3.velocity + k4.velocity);
      epoch.velocity    += (dt/6.) * (k1.acceleration + 2*k2.acceleration + 2*k3.acceleration + k4.acceleration);
      epoch.acceleration = acceleration(epoch, forces, satellite, earthRotation, ephemerides);
      orbit.push_back(epoch);
    }, timing);

    return orbit;
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif /* __GROOPS_ORBITPROPAGATORRUNGEKUTTA4__ */
