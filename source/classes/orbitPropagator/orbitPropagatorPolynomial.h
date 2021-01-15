/***********************************************/
/**
* @file orbitPropagatorPolynomial.h
*
* @brief Propagate a dynamic orbit using an integration Polynomial method.
* @see orbitPropagator
*
* @author Matthias Ellmer
* @date 2017-01-30
*
*/
/***********************************************/

#ifndef __GROOPS_ORBITPROPAGATORPOLYNOMIAL__
#define __GROOPS_ORBITPROPAGATORPOLYNOMIAL__

// Latex documentation
#ifdef DOCSTRING_OrbitPropagator
static const char *docstringOrbitPropagatorPolynomial = R"(
\subsection{Polynomial}\label{orbitPropagatorType:Polynomial}
This class implements an integration Polynomial method to propagate a satellite orbit under
the influence of \configClass{Forces}{forcesType}. Satellite is assumed to be oriented along-track.
Implementation is based on code by Torsten Mayer-Gürr.
)";
#endif

/***********************************************/

#include "classes/orbitPropagator/orbitPropagator.h"

/***** CLASS ***********************************/

/** @brief Propagate orbit using Polynomial method.
* @ingroup orbitPropagatorGroup
* @see orbitPropagator */
class OrbitPropagatorPolynomial : public OrbitPropagator
{
  OrbitPropagatorPtr warmup;   // Used for generation of warmup values.
  UInt   degree;               // of integration polynomial.
  Int    shift;
  Double epsilon;
  Matrix factorPos, factorVel;
  Bool   corrector;

public:
  OrbitPropagatorPolynomial(Config &config);

  OrbitArc integrateArc(OrbitEpoch startEpoch, Time sampling, UInt posCount, ForcesPtr forces, SatelliteModelPtr satellite,
                        EarthRotationPtr earthRotation, EphemeridesPtr ephemerides, Bool timing) const override;
};

/***********************************************/

inline OrbitPropagatorPolynomial::OrbitPropagatorPolynomial(Config &config)
{
  try
  {
    readConfig(config, "degree",    degree,    Config::MUSTSET, "9",           "polynomial degree to integrate accelerations");
    readConfig(config, "shift",     shift,     Config::MUSTSET, "2",           "shift polynomial in future (predicted accelerations)");
    readConfig(config, "epsilon",   epsilon,   Config::MUSTSET, "1e-6",        "[m] max. position change to recompute forces");
    readConfig(config, "warmup",    warmup,    Config::MUSTSET, "rungeKutta4", "to compute epochs before start epoch");
    readConfig(config, "corrector", corrector, Config::DEFAULT,  "0",          "apply corrector iteration if position change is larger than epsilon");
    if (isCreateSchema(config)) return;

    if(!(degree%2))
      throw(Exception("Polynomial Propagator: Polynomial degree must be odd. Is: ["+ degree % "%i"s +"]"));

    if(shift<0)
      throw(Exception("Polynomial Propagator: Polynomial shift must be positive. Is: ["+ shift % "%i"s +"]"));

    if(shift>static_cast<Int>(degree))
      throw(Exception("Polynomial Propagator: Polynomial shift must be smaller than degree. Is: ["+ shift % "%i"s +"]"));

    // Solve equation system for integration polynomial
    // tau = 0 at i = shift-degree
    Matrix  W(degree+1, degree+1);
    for(UInt i=0; i<=degree; i++)
      for(UInt n=0; n<=degree; n++)
        W(i,n) = ((n==0) ? 1. : pow(static_cast<Double>(i)+shift-degree, n));
    inverse(W);

    factorPos = factorVel = Matrix(shift+1, degree+1); // shift+1 epochs integration in future
    for(UInt i=0; i<factorPos.rows(); i++)
      for(UInt k=0; k<=degree; k++)
        for(UInt n=0; n<=degree; n++)
        {
          factorPos(i,k) += 1./((n+1)*(n+2)) * pow(i+1, n+2) * W(n,k);
          factorVel(i,k) += 1./(n+1)         * pow(i+1, n+1) * W(n,k);
        }
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline OrbitArc OrbitPropagatorPolynomial::integrateArc(OrbitEpoch startEpoch, Time sampling, UInt posCount, ForcesPtr forces,
                                                       SatelliteModelPtr satellite, EarthRotationPtr earthRotation, EphemeridesPtr ephemerides, Bool timing) const
{
  try
  {
    // Compute warmup
    OrbitArc orbit = flip(warmup->integrateArc(startEpoch, -sampling, degree-shift+1, forces, satellite, earthRotation, ephemerides, FALSE));
    orbit.append( warmup->integrateArc(startEpoch, sampling, shift+1, forces, satellite, earthRotation, ephemerides, FALSE).subArc(1,shift) );
    const UInt idx0 = degree-shift; // index of startEpoch

    // Actual integration
    const Double dt = sampling.seconds();
    UInt maxCorrectors = 10;

    Single::forEach(posCount-1, [&](UInt k)
    {
      // Integration: #shift epochs refinement, one new epoch prediction
      OrbitEpoch epoch;
      epoch.time = orbit.back().time + sampling;
      orbit.push_back(epoch);

      Bool correctorRepeatLoop = FALSE;
      UInt correctorIterations = 0;
      do
      {
        correctorRepeatLoop = FALSE;
        for(UInt i=0; i<factorPos.rows(); i++)
        {
          const Vector3d posOld = orbit.at(idx0+k+1+i).position;
          orbit.at(idx0+k+1+i).position = orbit.at(idx0+k).position + (i+1)*dt * orbit.at(idx0+k).velocity;
          orbit.at(idx0+k+1+i).velocity = orbit.at(idx0+k).velocity;
          for(UInt n=0; n<=degree; n++)
          {
            orbit.at(idx0+k+1+i).position += dt*dt*factorPos(i,n) * orbit.at(idx0+k+n+shift-degree).acceleration;
            orbit.at(idx0+k+1+i).velocity +=    dt*factorVel(i,n) * orbit.at(idx0+k+n+shift-degree).acceleration;
          }
          // recompute forces
          if((orbit.at(idx0+k+1+i).position-posOld).r()>epsilon)
          {
            orbit.at(idx0+k+1+i).acceleration = acceleration(orbit.at(idx0+k+1+i), forces, satellite, earthRotation, ephemerides);
            correctorRepeatLoop = corrector && TRUE;
            correctorIterations++;
          }
        }
      } while(correctorRepeatLoop && (correctorIterations < maxCorrectors));
    }, timing);

    // remove warmup and last epochs
    orbit.remove(0, degree-shift);
    orbit.remove(orbit.size()-shift, shift);

    return orbit;
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif /* __GROOPS_ORBITPROPAGATORPOLYNOMIAL__ */
