/***********************************************/
/**
* @file miscAccelerationsRelativisticEffect.h
*
* @brief Relativistic effect (IERS2010).
* @see MiscAccelerations
*
* @author Torsten Mayer-Guerr
* @date 2014-03-15
*
*/
/***********************************************/

#ifndef __GROOPS_MISCACCELERATIONSRELATIVISTICEFFECT__
#define __GROOPS_MISCACCELERATIONSRELATIVISTICEFFECT__

// Latex documentation
#ifdef DOCSTRING_MiscAccelerations
static const char *docstringMiscAccelerationsRelativisticEffect = R"(
\subsection{Relativistic effect}\label{miscAccelerationsType:relativisticEffect}
The relativistic effect to the acceleration of an artificial Earth satellite
according to IERS2010 conventions.
)";
#endif

/***********************************************/

#include "base/planets.h"
#include "inputOutput/logging.h"
#include "classes/miscAccelerations/miscAccelerations.h"

/***** CLASS ***********************************/

/** @brief Relativistic effect (IERS2010).
* @ingroup miscAccelerationsGroup
* @see MiscAccelerations */
class MiscAccelerationsRelativisticEffect : public MiscAccelerationsBase
{
  Double   factor;
  Double   GM;
  Double   beta, gamma;
  Vector3d J;

public:
  MiscAccelerationsRelativisticEffect(Config &config);

  Vector3d acceleration(SatelliteModelPtr satellite, const Time &time, const Vector3d &position, const Vector3d &velocity,
                        const Rotary3d &rotSat, const Rotary3d &rotEarth, EphemeridesPtr ephemerides) override;
};

/***********************************************/

inline MiscAccelerationsRelativisticEffect::MiscAccelerationsRelativisticEffect(Config &config)
{
  try
  {
    Double momentum;

    readConfig(config, "beta",   beta,     Config::DEFAULT,  "1.0",   "PPN (parameterized post-Newtonian) parameter");
    readConfig(config, "gamma",  gamma,    Config::DEFAULT,  "1.0",   "PPN (parameterized post-Newtonian) parameter");
    readConfig(config, "J",      momentum, Config::DEFAULT,  "9.8e8", "Earth’s angular momentum per unit mass [m**2/s]");
    readConfig(config, "GM",     GM,       Config::DEFAULT,  STRING_DEFAULT_GM, "Geocentric gravitational constant");
    readConfig(config, "factor", factor,   Config::DEFAULT,  "1.0",   "the result is multplied by this factor");
    if(isCreateSchema(config)) return;

    J = Vector3d(0, 0, momentum);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Vector3d MiscAccelerationsRelativisticEffect::acceleration(SatelliteModelPtr /*satellite*/, const Time &time,
                                                                  const Vector3d &pos, const Vector3d &vel,
                                                                  const Rotary3d &/*rotSat*/, const Rotary3d &rotEarth, EphemeridesPtr ephemerides)
{
  try
  {
    if(!ephemerides)
      throw(Exception("No ephemerides given"));
    if(vel.r() == 0.)
      logWarning<<"MiscAccelerationsRelativisticEffect: velocity is zero"<<Log::endl;

    // Position and velocity of sun relative to Earth
    Vector3d posSun, velSun;
    ephemerides->ephemeris(time, Ephemerides::SUN, posSun, velSun);

    const Double   r = pos.r();
    const Vector3d a = GM/pow(LIGHT_VELOCITY,2)/pow(r,3) * ((2*(beta+gamma)*GM/r-gamma*vel.quadsum())*pos + 2*(1+gamma)*inner(pos,vel) * vel
                                                           +(1+gamma)*(3/r/r*inner(pos,J)*crossProduct(pos,vel) + crossProduct(vel,J)))
                     + GM_Sun/pow(LIGHT_VELOCITY,2)/pow(posSun.r(),3) * (1+2*gamma) * crossProduct(crossProduct(velSun, posSun), vel);
    return factor * rotEarth.rotate(a);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
