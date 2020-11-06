/***********************************************/
/**
* @file tidesAstronomical.cpp
*
* @brief Astronomical tides
* Computed from planetary ephemerides.
* @see Tides
*
* @author Torsten Mayer-Guerr
* @date 2015-01-24
*
*/
/***********************************************/

#include "base/import.h"
#include "base/sphericalHarmonics.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/tides/tides.h"
#include "classes/tides/tidesAstronomical.h"

/***********************************************/

TidesAstronomical::TidesAstronomical(Config &config)
{
  try
  {
    Double c20;
    std::string choice;

    readConfig(config, "useMoon",    useMoon,    Config::DEFAULT,  "1", "TGP of moon");
    readConfig(config, "useSun",     useSun,     Config::DEFAULT,  "1", "TGP of sun");
    readConfig(config, "usePlanets", usePlanets, Config::DEFAULT,  "1", "TGP of planets");
    readConfig(config, "useEarth",   useEarth,   Config::DEFAULT,  "1", "TGP of Earth");
    readConfig(config, "c20Earth",   c20,        Config::DEFAULT,  "-4.84166854896119e-04", "J2 flattening of the Earth");
    readConfig(config, "factor",     factor,     Config::DEFAULT,  "1.0", "the result is multplied by this factor, set -1 to substract the field");
    if(isCreateSchema(config)) return;

    // Earth with flattening
    // ---------------------
    Matrix cnm(3,Matrix::TRIANGULAR, Matrix::LOWER);
    Matrix snm(3,Matrix::TRIANGULAR, Matrix::LOWER);
    cnm(2,0) = c20;
    j2earth = SphericalHarmonics(GM_Earth, R_Earth, cnm, snm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/************************************************/

// Tidal potential at point2 of a body with point mass GM at point1
// relative to coordinate center
Double TidesAstronomical::directTidePotential(Double GM, const Vector3d &point1, const Vector3d &point2) const
{
  Vector3d diff = point2-point1;
  Double r1  = point1.r();
  Double r12 = diff.r();

  return -GM*(1/r1 + 1/(r1*r1*r1)*inner(point1,point2) - 1/r12);
}

/************************************************/

Double TidesAstronomical::directTideRadialGradient(Double /*GM*/, const Vector3d &/*point1*/, const Vector3d &/*point2*/) const
{
  try
  {
    throw(Exception("not yet implemented"));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/************************************************/

// Tidal acceleration at point2 of a body with point mass GM at point1
// relative to coordinate center
Vector3d TidesAstronomical::directTideAcceleration(Double GM, const Vector3d &point1, const Vector3d &point2) const
{
  Vector3d diff = point2-point1;
  Double   r12  = diff.r();
  Double   r1   = point1.r();

  return -GM/(r12*r12*r12)*diff - GM/(r1*r1*r1)*point1;
}

/************************************************/

// Tidal acceleration gradient at point2 of a body with point mass GM at point1
// relative to coordinate center
Tensor3d TidesAstronomical::directTideGradient(Double GM, const Vector3d &point1, const Vector3d &point2) const
{
  Vector3d diff = point2-point1;
  Double   r12 = diff.r();
  Double   r3  = pow(r12,3);
  Double   r5  = pow(r12,5);
  Tensor3d T;
  T.xx() = 3*GM*diff.x()*diff.x()/r5 - 1/r3;
  T.yy() = 3*GM*diff.y()*diff.y()/r5 - 1/r3;
  T.zz() = 3*GM*diff.y()*diff.y()/r5 - 1/r3;
  T.xy() = 3*GM*diff.x()*diff.y()/r5;
  T.xz() = 3*GM*diff.x()*diff.z()/r5;
  T.yz() = 3*GM*diff.y()*diff.z()/r5;

  return T;
}

/************************************************/

Double TidesAstronomical::potential(const Time &time, const Vector3d &point, const Rotary3d &rotEarth, EarthRotationPtr /*rotation*/, EphemeridesPtr ephemerides) const
{
  try
  {
    if(!ephemerides)
      throw(Exception("No ephemerides given"));

    Vector3d posSat = rotEarth.inverseRotate(point);

    Double V = 0;
    if(useEarth && (ephemerides->origin() != Ephemerides::EARTH))
      V += directTidePotential(GM_Sun,  ephemerides->position(time, Ephemerides::EARTH), posSat);
    if(useMoon && (ephemerides->origin() != Ephemerides::MOON))
      V += directTidePotential(GM_Moon, ephemerides->position(time, Ephemerides::MOON),  posSat);
    if(useSun && (ephemerides->origin() != Ephemerides::SUN))
      V += directTidePotential(GM_Sun,  ephemerides->position(time, Ephemerides::SUN),   posSat);
    if(usePlanets)
    {
      V += directTidePotential(GM_MERCURY, ephemerides->position(time, Ephemerides::MERCURY), posSat)  // Mercury
         + directTidePotential(GM_VENUS  , ephemerides->position(time, Ephemerides::VENUS),   posSat)  // Venus
         + directTidePotential(GM_MARS   , ephemerides->position(time, Ephemerides::MARS),    posSat)  // Mars
         + directTidePotential(GM_JUPITER, ephemerides->position(time, Ephemerides::JUPITER), posSat)  // Jupiter
         + directTidePotential(GM_SATURN , ephemerides->position(time, Ephemerides::SATURN),  posSat); // Saturn
    }
    return factor*V;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/************************************************/

Double TidesAstronomical::radialGradient(const Time &time, const Vector3d &point, const Rotary3d &rotEarth, EarthRotationPtr /*rotation*/, EphemeridesPtr ephemerides) const
{
  try
  {
    if(!ephemerides)
      throw(Exception("No ephemerides given"));

    Vector3d posSat = rotEarth.inverseRotate(point);

    Double dVdr = 0;
    if(useEarth && (ephemerides->origin() != Ephemerides::EARTH))
      dVdr += directTideRadialGradient(GM_Sun,  ephemerides->position(time, Ephemerides::EARTH), posSat);
    if(useMoon && (ephemerides->origin() != Ephemerides::MOON))
      dVdr += directTideRadialGradient(GM_Moon, ephemerides->position(time, Ephemerides::MOON),  posSat);
    if(useSun && (ephemerides->origin() != Ephemerides::SUN))
      dVdr += directTideRadialGradient(GM_Sun,  ephemerides->position(time, Ephemerides::SUN),   posSat);
    if(usePlanets)
    {
      dVdr += directTideRadialGradient(GM_MERCURY, ephemerides->position(time, Ephemerides::MERCURY), posSat)  // Mercury
            + directTideRadialGradient(GM_VENUS  , ephemerides->position(time, Ephemerides::VENUS),   posSat)  // Venus
            + directTideRadialGradient(GM_MARS   , ephemerides->position(time, Ephemerides::MARS),    posSat)  // Mars
            + directTideRadialGradient(GM_JUPITER, ephemerides->position(time, Ephemerides::JUPITER), posSat)  // Jupiter
            + directTideRadialGradient(GM_SATURN , ephemerides->position(time, Ephemerides::SATURN),  posSat); // Saturn
    }
    return factor*dVdr;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector3d TidesAstronomical::gravity(const Time &time, const Vector3d &point, const Rotary3d &rotEarth, EarthRotationPtr /*rotation*/, EphemeridesPtr ephemerides) const
{
  try
  {
    if(!ephemerides)
      throw(Exception("No ephemerides given"));

    Vector3d posSat = rotEarth.inverseRotate(point);

    Vector3d g;
    if(useEarth && (ephemerides->origin() != Ephemerides::EARTH))
    {
      Vector3d posEarth = ephemerides->position(time, Ephemerides::EARTH);
      g += directTideAcceleration(GM_Earth, posEarth, posSat);
      g += j2earth.gravity(posSat-posEarth) - j2earth.gravity(-posEarth);
    }
    if(useMoon && (ephemerides->origin() != Ephemerides::MOON))
    {
      Vector3d posMoon = ephemerides->position(time, Ephemerides::MOON);
      g += directTideAcceleration(GM_Moon, posMoon, posSat);
      if(ephemerides->origin() == Ephemerides::EARTH)
        g -= GM_Moon/GM_Earth * j2earth.gravity(-posMoon);
    }
    if(useSun && (ephemerides->origin() != Ephemerides::SUN))
      g += directTideAcceleration(GM_Sun,  ephemerides->position(time, Ephemerides::SUN), posSat);
    if(usePlanets)
    {
      g += directTideAcceleration(GM_MERCURY, ephemerides->position(time, Ephemerides::MERCURY), posSat)  // Mercury
         + directTideAcceleration(GM_VENUS  , ephemerides->position(time, Ephemerides::VENUS),   posSat)  // Venus
         + directTideAcceleration(GM_MARS   , ephemerides->position(time, Ephemerides::MARS),    posSat)  // Mars
         + directTideAcceleration(GM_JUPITER, ephemerides->position(time, Ephemerides::JUPITER), posSat)  // Jupiter
         + directTideAcceleration(GM_SATURN , ephemerides->position(time, Ephemerides::SATURN),  posSat); // Saturn
    }
    return factor*rotEarth.rotate(g);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Tensor3d TidesAstronomical::gravityGradient(const Time &time, const Vector3d &point, const Rotary3d &rotEarth, EarthRotationPtr /*rotation*/, EphemeridesPtr ephemerides) const
{
  try
  {
    if(!ephemerides)
      throw(Exception("No ephemerides given"));

    Vector3d posSat = rotEarth.inverseRotate(point);

    Tensor3d T;
    if(useEarth && (ephemerides->origin() != Ephemerides::EARTH))
      T += directTideGradient(GM_Sun,  ephemerides->position(time, Ephemerides::EARTH), posSat);
    if(useMoon && (ephemerides->origin() != Ephemerides::MOON))
      T += directTideGradient(GM_Moon, ephemerides->position(time, Ephemerides::MOON),  posSat);
    if(useSun && (ephemerides->origin() != Ephemerides::SUN))
      T += directTideGradient(GM_Sun,  ephemerides->position(time, Ephemerides::SUN),   posSat);
    if(usePlanets)
    {
      T += directTideGradient(GM_MERCURY, ephemerides->position(time, Ephemerides::MERCURY), posSat)  // Mercury
         + directTideGradient(GM_VENUS  , ephemerides->position(time, Ephemerides::VENUS),   posSat)  // Venus
         + directTideGradient(GM_MARS   , ephemerides->position(time, Ephemerides::MARS),    posSat)  // Mars
         + directTideGradient(GM_JUPITER, ephemerides->position(time, Ephemerides::JUPITER), posSat)  // Jupiter
         + directTideGradient(GM_SATURN , ephemerides->position(time, Ephemerides::SATURN),  posSat); // Saturn
    }
    return factor*rotEarth.rotate(T);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector3d TidesAstronomical::deformation(const Time &/*time*/, const Vector3d &/*point*/, const Rotary3d &/*rotEarth*/,
                                        EarthRotationPtr /*rotation*/, EphemeridesPtr /*ephemerides*/,
                                        Double /*gravity*/, const Vector &/*hn*/, const Vector &/*ln*/) const
{
  return Vector3d(0,0,0);
}

/***********************************************/

void TidesAstronomical::deformation(const std::vector<Time> &/*time*/, const std::vector<Vector3d> &/*point*/, const std::vector<Rotary3d> &/*rotEarth*/,
                              EarthRotationPtr /*rotation*/, EphemeridesPtr /*ephemerides*/, const std::vector<Double> &/*gravity*/, const Vector &/*hn*/, const Vector &/*ln*/,
                              std::vector<std::vector<Vector3d>> &/*disp*/) const
{
}

/***********************************************/

SphericalHarmonics TidesAstronomical::sphericalHarmonics(const Time &/*time*/, const Rotary3d &/*rotEarth*/, EarthRotationPtr /*rotation*/, EphemeridesPtr /*ephemerides*/,
                                                         UInt /*maxDegree*/, UInt /*minDegree*/, Double /*GM*/, Double /*R*/) const
{
  // tide generating potential is not harmonic
  return SphericalHarmonics();
}

/***********************************************/
/***********************************************/
