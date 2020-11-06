/***********************************************/
/**
* @file eclipseSOLAARS.h
*
* @brief Shadowing of satellites by moon and Earth.
* @see Eclipse
*
* @author Torsten Mayer-Guerr
* @date 2020-03-08
*
*/
/***********************************************/

#ifndef __GROOPS_ECLIPSESOLAARS__
#define __GROOPS_ECLIPSESOLAARS__

// Latex documentation
#ifdef DOCSTRING_Eclipse
static const char *docstringEclipseSOLAARS = R"(
\subsection{SOLAARS}
Earth’s penumbra modeling with Solar radiation pressure with
Oblateness and Lower Atmospheric Absorption, Refraction, and Scattering (SOLAARS).
See Robertson, Robbie. (2015),
Highly Physical Solar Radiation Pressure Modeling During Penumbra Transitions (pp. 67-75).
)";
#endif

/***********************************************/

#include "classes/eclipse/eclipse.h"

/***** CLASS ***********************************/

/** @brief Shadowing of satellites by moon and Earth.
* @ingroup eclipseGroup
* @see Eclipse */
class EclipseSOLAARS : public Eclipse
{
public:
  EclipseSOLAARS(Config &/*config*/) {}

  virtual Double factor(const Time &timeGPS, const Vector3d &position, EphemeridesPtr ephemerides) const override;
};

/***********************************************/

inline Double EclipseSOLAARS::factor(const Time &timeGPS, const Vector3d &position, EphemeridesPtr ephemerides) const
{
  try
  {
    if(!ephemerides)
      throw(Exception("No ephemerides given"));

    const Vector3d posSun = ephemerides->position(timeGPS, Ephemerides::SUN);
    Double factorSun      = shadowScalingFactor(position, posSun, Vector3d(), R_Earth+100e3);
    Double factorMoon     = shadowScalingFactor(position, posSun, ephemerides->position(timeGPS, Ephemerides::MOON), R_Moon);
    if(factorSun == 1.)
      return factorSun * factorMoon;

    // For detailed information on this calculation, see Robertson (2015, pp. 67-75)
    const Vector3d R  = normalize(posSun);
    const Double   rR = -1e-06 * inner(position, R); // in 1e6 meters
    Vector3d r = position;
    r.z() *= 6378137./6356752.; // account for flattening
    const Double rE = 1e-6 * (r-R*inner(r, R)).r();

    // model parameters
    constexpr Double b11 =  0.1715;   constexpr Double b12 = -0.1423;    constexpr Double b13 =  0.01061; constexpr Double b14 = -0.01443;
    constexpr Double b21 =  0.008162; constexpr Double b22 =  0.3401;
    constexpr Double b31 =  260.9;    constexpr Double b32 = -0.4661;    constexpr Double b33 =  27.81;   constexpr Double b34 = -0.009437;
    constexpr Double b41 = -0.006119; constexpr Double b42 =  1.176;     constexpr Double b43 =  6.385;
    constexpr Double b51 =  87.56;    constexpr Double b52 = -0.09188;   constexpr Double b53 =  19.30;   constexpr Double b54 = -0.01089;
    constexpr Double b61 =  0.002047; constexpr Double b62 =  6.409;
    constexpr Double b71 =  61.98;    constexpr Double b72 = -0.1629;    constexpr Double b73 =  27.87;   constexpr Double b74 = -0.02217;
    constexpr Double b81 =  6.413;    constexpr Double b82 = -0.0002593; constexpr Double b83 = -0.01479; constexpr Double b84 = -0.1318;

    const Double a1 = b11 * std::exp(b12 * rR) + b13 * std::exp(b14 * rR);
    const Double a2 = b21 * rR                 + b22;
    const Double a3 = b31 * std::exp(b32 * rR) + b33 * std::exp(b34 * rR);
    const Double a4 = b41 * std::pow(rR, b42)  + b43;
    const Double a5 = b51 * std::exp(b52 * rR) + b53 * std::exp(b54 * rR);
    const Double a6 = b61 * rR                 + b62;
    const Double a7 = b71 * std::exp(b72 * rR) + b73 * std::exp(b74 * rR);
    const Double a8 = b81 * std::exp(b82 * rR) + b83 * std::exp(b84 * rR);

    factorSun = 0.5*(1+a1+a2+a1*std::tanh(a3*(rE-a4))+a2*std::tanh(a5*(rE-a6))+std::tanh(a7*(rE-a8)))/(1+a1+a2);
    if(std::isnan(factorSun))
      factorSun = 1;
    return factorSun * factorMoon;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
