/***********************************************/
/**
* @file orbitPropagator.cpp
*
* @brief Propagate orbit using different algorithms.
*
* @author Matthias Ellmer
* @date 2017-01-19
*
*/
/***********************************************/

#define DOCSTRING_OrbitPropagator

#include "base/import.h"
#include "config/configRegister.h"
#include "inputOutput/logging.h"
#include "classes/orbitPropagator/orbitPropagatorEuler.h"
#include "classes/orbitPropagator/orbitPropagatorRungeKutta4.h"
#include "classes/orbitPropagator/orbitPropagatorAdamsBashforthMoulton.h"
#include "classes/orbitPropagator/orbitPropagatorStoermerCowell.h"
#include "classes/orbitPropagator/orbitPropagatorGaussJackson.h"
#include "classes/orbitPropagator/orbitPropagatorPolynomial.h"
#include "classes/orbitPropagator/orbitPropagatorFile.h"
#include "classes/orbitPropagator/orbitPropagator.h"

/***********************************************/

GROOPS_REGISTER_CLASS(OrbitPropagator, "orbitPropagatorType",
                      OrbitPropagatorEuler,
                      OrbitPropagatorRungeKutta4,
                      OrbitPropagatorAdamsBashforthMoulton,
                      OrbitPropagatorStoermerCowell,
                      OrbitPropagatorGaussJackson,
                      OrbitPropagatorPolynomial,
                      OrbitPropagatorFile)

GROOPS_READCONFIG_CLASS(OrbitPropagator, "orbitPropagatorType")

/***********************************************/

OrbitPropagatorPtr OrbitPropagator::create(Config &config, const std::string &name)
{
  try
  {
    OrbitPropagatorPtr orbitPropagator;
    std::string choice;

    readConfigChoice(config, name, choice, Config::MUSTSET, "", "algorithm for integration");
    if (readConfigChoiceElement(config, "euler",               choice, "Euler step method"))
      orbitPropagator = OrbitPropagatorPtr(new OrbitPropagatorEuler(config));
    if (readConfigChoiceElement(config, "rungeKutta4",         choice, "The classical Runge-Kutta 4 method"))
      orbitPropagator = OrbitPropagatorPtr(new OrbitPropagatorRungeKutta4(config));
    if (readConfigChoiceElement(config, "adamsBashforthMoulton", choice, "Adams-Bashforth-Moulton class of predictor-corrector method"))
      orbitPropagator = OrbitPropagatorPtr(new OrbitPropagatorAdamsBashforthMoulton(config));
    if (readConfigChoiceElement(config, "stoermerCowell",      choice, "Stoermer-Cowell predictor-corrector method"))
      orbitPropagator = OrbitPropagatorPtr(new OrbitPropagatorStoermerCowell(config));
    if (readConfigChoiceElement(config, "gaussJackson",        choice, "Gauss-Jackson method"))
      orbitPropagator = OrbitPropagatorPtr(new OrbitPropagatorGaussJackson(config));
    if (readConfigChoiceElement(config, "polynomial",          choice, "Integration Polynomial method using prediction only"))
      orbitPropagator = OrbitPropagatorPtr(new OrbitPropagatorPolynomial(config));
    if (readConfigChoiceElement(config, "file",                choice, "Read orbit from file"))
      orbitPropagator = OrbitPropagatorPtr(new OrbitPropagatorFile(config));
    endChoice(config);

    return orbitPropagator;
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Rotary3d OrbitPropagator::orientation(const Time &/*time*/, const Vector3d &position, const Vector3d &velocity, SatelliteModelPtr /*satellite*/) const
{
  return Rotary3d(velocity, crossProduct(velocity, position));
}

/***********************************************/

Arc OrbitPropagator::flip(const Arc &arc)
{
  try
  {
    Arc arcNew;
    for(UInt i=arc.size(); i-->0;)
      arcNew.push_back(arc.at(i));
    return arcNew;
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector3d OrbitPropagator::acceleration(const OrbitEpoch &epoch, ForcesPtr forces, SatelliteModelPtr satellite,
                                       EarthRotationPtr earthRotation, EphemeridesPtr ephemerides) const
{
  try
  {
    const Rotary3d rotSat   = orientation(epoch.time, epoch.position, epoch.velocity, satellite);
    const Rotary3d rotEarth = earthRotation->rotaryMatrix(epoch.time);
    return rotEarth.inverseRotate(forces->acceleration(satellite, epoch.time, epoch.position, epoch.velocity, rotSat, rotEarth, earthRotation, ephemerides));
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
