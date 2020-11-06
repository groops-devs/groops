/***********************************************/
/**
* @file orbitPropagatorEuler.h
*
* @brief Propagate a dynamic orbit using Euler's method.
* @see orbitPropagator
*
* @author Matthias Ellmer
* @date 2017-01-19
*
*/
/***********************************************/

#ifndef __GROOPS_ORBITPROPAGATOREULER__
#define __GROOPS_ORBITPROPAGATOREULER__

// Latex documentation
#ifdef DOCSTRING_OrbitPropagator
static const char *docstringOrbitPropagatorEuler = R"(
\subsection{Euler}
This class implements Euler's method to propagate a satellite orbit under the influence of \configClass{Forces}{forcesType}.
Satellite is assumed to be oriented along-track.
)";
#endif

/***********************************************/

#include "classes/orbitPropagator/orbitPropagator.h"

/***** CLASS ***********************************/

/** @brief Propagate orbit using Euler's method.
* @ingroup orbitPropagatorGroup
* @see orbitPropagator */
class OrbitPropagatorEuler : public OrbitPropagator
{
public:
  OrbitPropagatorEuler(Config &/*config*/) {}

  OrbitArc integrateArc(OrbitEpoch startEpoch, Time sampling, UInt posCount, ForcesPtr forces, SatelliteModelPtr satellite, EarthRotationPtr earthRotation,
                        EphemeridesPtr ephemerides, Bool timing) const override;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline OrbitArc OrbitPropagatorEuler::integrateArc(OrbitEpoch startEpoch, Time sampling, UInt posCount, ForcesPtr forces,
                                                   SatelliteModelPtr satellite, EarthRotationPtr earthRotation, EphemeridesPtr ephemerides, Bool timing) const
{
  try
  {
    OrbitArc orbit;
    startEpoch.acceleration = acceleration(startEpoch, forces, satellite, earthRotation, ephemerides);
    orbit.push_back(startEpoch);
    const Double dt = sampling.seconds();

    if(timing) logTimerStart;
    for(UInt k=1; k<posCount; k++)
    {
      if(timing) logTimerLoop(k, posCount);
      OrbitEpoch epoch;
      epoch.time         = orbit.at(k-1).time + sampling;
      epoch.position     = orbit.at(k-1).position + dt * orbit.at(k-1).velocity;
      epoch.velocity     = orbit.at(k-1).velocity + dt * orbit.at(k-1).acceleration;
      epoch.acceleration = acceleration(epoch, forces, satellite, earthRotation, ephemerides);
      orbit.push_back(epoch);
    }
    if(timing) logTimerLoopEnd(posCount);

    return orbit;
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif /* __GROOPS_ORBITPROPAGATOREULER__ */
