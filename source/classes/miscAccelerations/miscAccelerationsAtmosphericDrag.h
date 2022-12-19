/***********************************************/
/**
* @file miscAccelerationsAtmosphericDrag.h
*
* @brief Atmospheric drag.
* @see MiscAccelerations
*
* @author Beate Klinger
* @author Sandro Krauss
* @date 2015-12-10
*
*/
/***********************************************/

#ifndef __GROOPS_MISCACCELERATIONSATMOSPHERICDRAG__
#define __GROOPS_MISCACCELERATIONSATMOSPHERICDRAG__

// Latex documentation
#ifdef DOCSTRING_MiscAccelerations
static const char *docstringMiscAccelerationsAtmosphericDrag = R"(
\subsection{AtmosphericDrag}\label{miscAccelerationsType:atmosphericDrag}
Atmospheric drag model.
Algorithm for the atmospheric drag modelling is based on the free molecule flow
theory by Sentman 1961. An analytical expression of this treatise is given in
Moe and Moe 2005.

Sentman L. (1961), Free molecule flow theory and its application to the determination
of aerodynamic forces, Technical report.

Moe K., Moe M. M. (2005), Gas-surface interactions and satellite drag coefficients,
Planetary and Space Science 53(8), 793-801, doi:10.1016/j.pss.2005.03.005.

Optional determination steps:
Turn temperature on or off.
In the first case, the model mentioned above is applied, which estimates variable drag
and lift coefficients - in the latter case a constant drag coefficient can be specified.

Turn wind on/off:
It enables the usage of the Horizontal Wind Model 2014 to add additional thermospheric
winds in the calculation process.
)";
#endif

/***********************************************/

#include "classes/thermosphere/thermosphere.h"
#include "classes/miscAccelerations/miscAccelerations.h"

/***** CLASS ***********************************/

/** @brief Atmospheric drag.
 * @ingroup miscAccelerationsGroup
 * @see MiscAccelerations */
class MiscAccelerationsAtmosphericDrag : public MiscAccelerationsBase
{
  ThermospherePtr thermosphere;
  Bool            useTemperature, useWind;
  Vector3d        omega;
  Double          factor;

public:
  MiscAccelerationsAtmosphericDrag(Config &config);

  Vector3d acceleration(SatelliteModelPtr satellite, const Time &time, const Vector3d &position, const Vector3d &velocity,
                        const Rotary3d &rotSat, const Rotary3d &rotEarth, EphemeridesPtr ephemerides) override;

  static Vector3d force(SatelliteModelPtr satellite, const Vector3d &direction, Double velocity, Double density, Double temperature);
};

/***********************************************/

inline MiscAccelerationsAtmosphericDrag::MiscAccelerationsAtmosphericDrag(Config &config)
{
  try
  {
    Double angleVelocity;

    readConfig(config, "thermosphere",        thermosphere,   Config::MUSTSET, "", "");
    readConfig(config, "earthRotation",       angleVelocity,  Config::DEFAULT, "7.29211585531e-5", "[rad/s]");
    readConfig(config, "considerTemperature", useTemperature, Config::DEFAULT, "1",   "compute drag and lift, otherwise simple drag coefficient is used");
    readConfig(config, "considerWind",        useWind,        Config::DEFAULT, "1",   "");
    readConfig(config, "factor",              factor,         Config::DEFAULT, "1.0", "the result is multiplied by this factor");
    if(isCreateSchema(config)) return;

    omega = Vector3d(0, 0, angleVelocity);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Vector3d MiscAccelerationsAtmosphericDrag::acceleration(SatelliteModelPtr satellite, const Time &time,
                                                               const Vector3d &position, const Vector3d &velocity,
                                                               const Rotary3d &rotSat, const Rotary3d &rotEarth, EphemeridesPtr /*ephemerides*/)
{
  try
  {
    if(!satellite)
      throw(Exception("No satellite model given"));

    Double   density, temperature;
    Vector3d wind;
    thermosphere->state(time, rotEarth.rotate(position), density, temperature, wind);
    if(!useTemperature)
      temperature = 0;
    if(!useWind)
      wind = Vector3d();

    // direction and speed of thermosphere relative to satellite in SRF
    Vector3d direction = rotSat.inverseRotate(rotEarth.inverseRotate(wind) + crossProduct(omega, position) - velocity);
    const Double v = direction.normalize();

    return (factor/satellite->mass) * rotEarth.rotate(rotSat.rotate(force(satellite, direction, v, density, temperature)));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Vector3d MiscAccelerationsAtmosphericDrag::force(SatelliteModelPtr satellite, const Vector3d &direction, Double velocity, Double density, Double temperature)
{
  try
  {
    constexpr Double R     = 8.3144;    // universal gas constant [J mol-1 K-1]
    constexpr Double M     = 15.5/1000; // average modelcular mass [kg/mol]
    constexpr Double alpha = 0.9;       // accomodation coefficientDrag [-]
    constexpr Double Tw    = 273;       // temperature of the surface [K]

    Vector3d F;
    for(UInt i=0; i<satellite->surfaces.size(); i++)
    {
      auto &surface = satellite->surfaces.at(i);
      if((surface.type == SatelliteModel::Surface::PLATE) || (surface.type == SatelliteModel::Surface::GRACESHADOW))
      {
        const Double cosPhi = -inner(direction, surface.normal);

        if(temperature > 0)
        {
          const Vector3d ul   = crossProduct(direction, crossProduct(direction, surface.normal));
          const Double   s    = velocity/std::sqrt(2*R*temperature/M); // molecular speed ratio
          const Double   P    = std::exp(-std::pow(cosPhi*s, 2))/s;
          const Double   G    = 1/(2*s*s);
          const Double   Z    = 1 + std::erf(cosPhi*s);
          const Double   x    = std::sqrt(1./6.*(1 + alpha * (3*R*Tw/(velocity*velocity)-1))) * (cosPhi*std::sqrt(PI)*Z+P);
          F += 0.5 * density * surface.area *velocity*velocity * ((P/std::sqrt(PI) + cosPhi*((1+G)*Z+x)) * direction + (G*Z+x) * ul);
        }
        else if(cosPhi > 0) // simple drag coefficient
          F += 0.5 * density * surface.area * cosPhi * velocity*velocity * satellite->coefficientDrag * direction;
      }
      else
        throw(Exception("sphere or cylinder not implemented yet"));
    } // for(i=surface)

    return F;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
