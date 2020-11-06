/***********************************************/
/**
* @file orbitPropagatorStoermerCowell.h
*
* @brief Propagate a dynamic orbit using an arbitrary order Stoermer-Cowell predictor-corrector algorithm.
* @see orbitPropagator
*
* @author Matthias Ellmer
* @date 2017-02-30
*
*/
/***********************************************/

#ifndef __GROOPS_ORBITPROPAGATORSTOERMERCOWELL__
#define __GROOPS_ORBITPROPAGATORSTOERMERCOWELL__

// Latex documentation
#ifdef DOCSTRING_OrbitPropagator
static const char *docstringOrbitPropagatorStoermerCowell = R"(
\subsection{StoermerCowell}
This class implements the Stoermer-Cowell class of predictor-corrector orbit propagators for a satellite orbit
under the influence of \configClass{Forces}{forcesType}. The coefficients for the Stoermer predictor and Cowell corrector
are derived using the equations given in section 4.2.6 of [1]. Stoermer-Cowell is a double integration algorithm,
yielding positions directly from accelertions. It does not produce velocities. The velocities are derived using
Adams-type propagators as suggested in [2]. Satellite is assumed to be oriented along-track.
[1] Montenbruck, Oliver, and Eberhard Gill. 2000. Satellite Orbits
[2] Berry, Matthew M., and Liam M. Healy. 2004. “Implementation of Gauss-Jackson Integration for Orbit Propagation.”
)";
#endif

/***********************************************/

#include "classes/orbitPropagator/orbitPropagator.h"

/***** CLASS ***********************************/

/** @brief Propagate orbit using Adams-Bashforth type method.
* @ingroup orbitPropagatorGroup
* @see orbitPropagator */
class OrbitPropagatorStoermerCowell : public OrbitPropagator
{
private:
  OrbitPropagatorPtr warmup;
  UInt   order;
  Vector stoermer, bashforth;
  Vector cowell,   moulton;

public:
  OrbitPropagatorStoermerCowell(Config &config);

  OrbitArc integrateArc(OrbitEpoch startEpoch, Time sampling, UInt posCount, ForcesPtr forces, SatelliteModelPtr satellite,
                        EarthRotationPtr earthRotation, EphemeridesPtr ephemerides, Bool timing) const override;

  /** @brief Compute Stoermer predictor coefficients.
  * @param order of the predictor
  * @returns @a Vector of length oder+1 */
  static Vector stoermerCoefficients(UInt order);

  /** @brief Compute Cowell corrector coefficients.
  * @param order of the corrector
  * @returns @a Vector of length oder+1 */
  static Vector cowellCoefficients(UInt order);

  /** @brief Compute backward difference operator coefficients.
  * @param degree of the backwards difference operator
  * @param index of the value to which the difference operator is applied
  * @param order of the integrator
  * @returns @a Vector of length @a order+1, with the backward difference factors
  * for the elements i = [0 .. order] of @a degree at the @a i-th position.
  *
  * backwardsDifference(0,1,2) -> [ 0.,  1.,  0.]
  * backwardsDifference(0,2,2) -> [ 0.,  0.,  1.]
  * backwardsDifference(1,2,2) -> [ 0., -1.,  1.]
  * backwardsDifference(2,2,2) -> [ 1., -2.,  1.]  */
  static Vector backwardsDifference(UInt degree, UInt index, UInt order);
};

/***********************************************/

inline OrbitPropagatorStoermerCowell::OrbitPropagatorStoermerCowell(Config &config)
{
  try
  {
    readConfig(config, "order",  order,  Config::MUSTSET, "4", "Order of the Stoermer-Cowell type propagator.");
    readConfig(config, "warmup", warmup, Config::MUSTSET, "rungeKutta4", "");
    if (isCreateSchema(config)) return;

    if(order<2)
      throw(Exception("Stoermer-Cowell must have an order of at least 2, is: "+order%"%i"s));

    stoermer = Vector(order+1); // Predictor for position
    const Vector coeffStoermer = stoermerCoefficients(order);
    for(UInt j=0; j<order; j++)
      axpy(coeffStoermer(j), backwardsDifference(j, order, order), stoermer);

    cowell = Vector(order+1);
    const Vector coeffCowell = cowellCoefficients(order);   // Corrector for position
    for(UInt j=0; j<order; j++)
      axpy(coeffCowell(j), backwardsDifference(j, order, order), cowell);

    // Set up coefficients for velocity
    bashforth = OrbitPropagatorAdamsBashforthMoulton::coefficients(OrbitPropagatorAdamsBashforthMoulton::factorsBashforth(order)); // Predictor for velocity
    moulton   = OrbitPropagatorAdamsBashforthMoulton::coefficients(OrbitPropagatorAdamsBashforthMoulton::factorsMoulton(order));   // Corrector for velocity
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline OrbitArc OrbitPropagatorStoermerCowell::integrateArc(OrbitEpoch startEpoch, Time sampling, UInt posCount, ForcesPtr forces,
                                                            SatelliteModelPtr satellite, EarthRotationPtr earthRotation, EphemeridesPtr ephemerides, Bool timing) const
{
  try
  {
    const Double dt = sampling.seconds();

    // Compute warmup
    OrbitArc orbit = flip(warmup->integrateArc(startEpoch, -sampling, order, forces, satellite, earthRotation, ephemerides, FALSE));

    // Integrate remaining arc
    if(timing) logTimerStart;
    for (UInt k=1; k<posCount; k++)
    {
      if(timing) logTimerLoop(k, posCount);

      // Predictor step
      // --------------
      OrbitEpoch epoch;
      epoch.time     = orbit.at(k+order-2).time + sampling;
      epoch.position += 2*orbit.at(k+order-2).position - orbit.at(k+order-3).position;
      for(UInt j=1; j<=order; j++) // Stoermer predictor
        epoch.position += dt*dt * stoermer(j) * orbit.at(k+j-2).acceleration;
      epoch.velocity += orbit.at(k+order-2).velocity;  // Adams-Bashforth for velocity component
      for(UInt j=1; j<=order; j++)
        epoch.velocity += dt * bashforth(j) * orbit.at(k+j-2).acceleration;
      epoch.acceleration = acceleration(epoch, forces, satellite, earthRotation, ephemerides);
      orbit.push_back(epoch);

      // Corrector step
      // --------------
      epoch.position = 2*orbit.at(k+order-2).position - orbit.at(k+order-3).position;
      for(UInt j=1; j<=order; j++) // Stoermer predictor
        epoch.position += dt*dt * cowell(j) * orbit.at(k+j-1).acceleration;
      epoch.velocity = orbit.at(k+order-2).velocity;  // Moulton corrector
      for(UInt j=1; j<=order; j++)
        epoch.velocity += dt * moulton(j) * orbit.at(k+j-1).acceleration;
      epoch.acceleration = acceleration(epoch, forces, satellite, earthRotation, ephemerides);
      orbit.at(k+order-1) = epoch;
    }
    if(timing) logTimerLoopEnd(posCount);

    // remove warmup
    orbit.remove(0, order-1);

    return orbit;
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

// Set up coefficients for explicit Stoermer predictor coefficients of requested order.
// Montenbruck, Oliver, and Eberhard Gill. 2000. Satellite Orbits, eq. 4.90
inline Vector OrbitPropagatorStoermerCowell::stoermerCoefficients(UInt order)
{
  Vector gamma = OrbitPropagatorAdamsBashforthMoulton::factorsMoulton(order);
  Vector stoermer(gamma.size());
  for(UInt j=0; j<gamma.size(); j++)
    stoermer(j) = (1.-j) * gamma(j);
  return stoermer;
}

/***********************************************/

// Set up coefficients for implicit Cowell corrector coefficients of requested order.
// Montenbruck, Oliver, and Eberhard Gill. 2000. Satellite Orbits, eq. 4.92
inline Vector OrbitPropagatorStoermerCowell::cowellCoefficients(UInt order)
{
  Vector stoermer = stoermerCoefficients(order);
  Vector cowell(order+1);
  cowell(0) = 1;
  for (UInt j = 1; j < cowell.size(); j++)
    cowell(j) = stoermer(j) - stoermer(j-1);
  return cowell;
}

/***********************************************/

// Montenbruck, Oliver, and Eberhard Gill. 2000. Satellite Orbits, eq. 4.53
inline Vector OrbitPropagatorStoermerCowell::backwardsDifference(UInt degree, UInt index, UInt order)
{
  if(degree == 0)
  {
    Vector a(order+1);
    a(index) = 1; // f(index) is 1 at the index-th position
    return a;
  }
  else if(degree == 1)
  {
    Vector a(order+1);
    a(index)   =  1; //  f(index)
    a(index-1) = -1; // -f(index-1)
    return a;
  }

  return backwardsDifference(degree-1, index, order) - backwardsDifference(degree-1, index-1, order);
}

/***********************************************/

#endif /* __GROOPS_ORBITPROPAGATORSTOERMERCOWELL__ */
