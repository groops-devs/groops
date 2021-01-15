/***********************************************/
/**
* @file orbitPropagatorAdamsBashforthMoulton.h
*
* @brief Propagate a dynamic orbit using Adams/Bashforth/Moulton algorithm.
* @see orbitPropagator
*
* @author Matthias Ellmer
* @date 2017-01-25
*
*/
/***********************************************/

#ifndef __GROOPS_ORBITPROPAGATORADAMSBASHFORTHMOULTON__
#define __GROOPS_ORBITPROPAGATORADAMSBASHFORTHMOULTON__

// Latex documentation
#ifdef DOCSTRING_OrbitPropagator
static const char *docstringOrbitPropagatorAdamsBashforthMoulton = R"(
\subsection{AdamsBashforthMoulton}
This class implements the Adams-Moulton class of predictor-corrector orbit propagators
for a satellite orbit under the influence of \configClass{Forces}{forcesType} using an implicit
Adams-Bashforth corrector. The coefficients for the propagator are derived using the equations
given in section 4.2.3 of [1]. Satellite is assumed to be oriented along-track.
[1] Montenbruck, Oliver, and Eberhard Gill. 2000. Satellite Orbits
)";
#endif

/***********************************************/

#include "classes/orbitPropagator/orbitPropagator.h"

/***** CLASS ***********************************/

/** @brief Propagate a dynamic orbit using Adams/Bashforth/Moulton algorithm.
* @ingroup orbitPropagatorGroup
* @see orbitPropagator */
class OrbitPropagatorAdamsBashforthMoulton : public OrbitPropagator
{
  OrbitPropagatorPtr warmup; // Used for generation of warmup values.
  Bool   applyMoultonCorrector;
  UInt   order;
  Vector betaAM, betaAB;

public:
  OrbitPropagatorAdamsBashforthMoulton(Config &config);

  OrbitArc integrateArc(OrbitEpoch startEpoch, Time sampling, UInt posCount, ForcesPtr forces, SatelliteModelPtr satellite,
                        EarthRotationPtr earthRotation, EphemeridesPtr ephemerides, Bool timing) const override;

  static Vector factorsBashforth(UInt order);
  static Vector factorsMoulton(UInt order);
  static Vector coefficients(const Vector &gamma);
};

/***********************************************/

inline OrbitPropagatorAdamsBashforthMoulton::OrbitPropagatorAdamsBashforthMoulton(Config &config)
{
  try
  {
    readConfig(config, "order",                  order,                  Config::MUSTSET,  "4", "Order of the Adams-Bashforth type propagator.");
    readConfig(config, "applyMoultonCorrector",  applyMoultonCorrector,  Config::DEFAULT,  "1", "Corrector step after Adams-Bashforth predcition.");
    readConfig(config, "warmup",                 warmup,                 Config::MUSTSET,  "rungeKutta4", "");
    if(isCreateSchema(config)) return;

    betaAB = coefficients(factorsBashforth(order)); // Adams-Bashforth predictor coefficients
    betaAM = coefficients(factorsMoulton(order));   // Coefficients for Adams-Moulton corrector
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline OrbitArc OrbitPropagatorAdamsBashforthMoulton::integrateArc(OrbitEpoch startEpoch, Time sampling, UInt posCount, ForcesPtr forces,
                                                                   SatelliteModelPtr satellite, EarthRotationPtr earthRotation, EphemeridesPtr ephemerides, Bool timing) const
{
  try
  {
    const Double dt = sampling.seconds();

    // Compute warmup
    OrbitArc orbit = flip(warmup->integrateArc(startEpoch, -sampling, order, forces, satellite, earthRotation, ephemerides, FALSE));

    // Integrate remaining arc
    Single::forEach(posCount-1, [&](UInt k)
    {
      // Predict using Adams-Bashforth
      OrbitEpoch epoch = orbit.back();
      epoch.time += sampling;
      for(UInt j=1; j<=order; j++)
        epoch.position += dt * betaAB(j) * orbit.at(k+j-1).velocity;
      for(UInt j=1; j<=order; j++)
        epoch.velocity += dt * betaAB(j) * orbit.at(k+j-1).acceleration;
      epoch.acceleration = acceleration(epoch, forces, satellite, earthRotation, ephemerides);
      orbit.push_back(epoch);

      // Correct using Adams-Moulton
      if(applyMoultonCorrector)
      {
        epoch = orbit.at(orbit.size()-2);
        epoch.time += sampling;
        for(UInt j=1; j<=order; j++)
          epoch.position += dt * betaAM(j) * orbit.at(k+j).velocity;
        for(UInt j=1; j<=order; j++)
          epoch.velocity += dt * betaAM(j) * orbit.at(k+j).acceleration;
        epoch.acceleration = acceleration(epoch, forces, satellite, earthRotation, ephemerides);
        orbit.at(k+order) = epoch;
      }
    }, timing);

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

// Set up coefficients for Adams-Bashforth predictor of requested order.
// Montenbruck, Oliver, and Eberhard Gill. 2000. Satellite Orbits, eq. 4.56
inline Vector OrbitPropagatorAdamsBashforthMoulton::factorsBashforth(UInt order)
{
  Vector gamma(order+1);
  for (UInt j=0; j<gamma.size(); j++)
  {
    gamma(j) = 1;
    for (UInt k=0; k<j; k++)
      gamma(j) -= (1./(1.+j-k)) * gamma(k);
  }
  return gamma;
}

/***********************************************/

// Set up coefficients for Adams-Moulton corrector of requested order.
// Montenbruck, Oliver, and Eberhard Gill. 2000. Satellite Orbits, eq. 4.64
inline Vector OrbitPropagatorAdamsBashforthMoulton::factorsMoulton(UInt order)
{
  Vector gamma(order+1);
  gamma(0) = 1;
  for (UInt j=1; j<gamma.size(); j++)
    for (UInt k=0; k<j; k++)
      gamma(j) -= (1./(1.+j-k)) * gamma(k);
  return gamma;
}

/***********************************************/

// Set up coefficients for Adams-Bashforth predictor of requested order.
// Montenbruck, Oliver, and Eberhard Gill. 2000. Satellite Orbits, eq. 4.60
inline Vector OrbitPropagatorAdamsBashforthMoulton::coefficients(const Vector &gamma)
{
  const UInt order = gamma.rows()-1;

  Matrix binomial(order, order, 1.);
  for(UInt i=2; i<binomial.rows(); i++)
    for(UInt k=1; k<i; k++)
      binomial(i,k) = binomial(i-1,k-1) + binomial(i-1,k);

  Vector beta(order+1);
  for(UInt j=1; j<beta.size(); j++)
    for(UInt l=order-j; l<order; l++)
      beta(j) += std::pow(-1, order-j) * gamma(l) * binomial(l, order-j);
  return beta;
}

/***********************************************/

#endif
