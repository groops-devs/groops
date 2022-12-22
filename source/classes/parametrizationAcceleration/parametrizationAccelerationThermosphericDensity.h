/***********************************************/
/**
* @file parametrizationAccelerationThermosphericDensity.h
*
* @brief Estimate the thermospheric density along the orbit.
*
* @author Torsten Mayer-Guerr
* @date 2020-11-06
*
*/
/***********************************************/

#ifndef __GROOPS_PARAMETRIZATIONTHERMOSPHERICDENSITY__
#define __GROOPS_PARAMETRIZATIONTHERMOSPHERICDENSITY__

// Latex documentation
#ifdef DOCSTRING_ParametrizationAcceleration
static const char *docstringParametrizationAccelerationThermosphericDensity = R"(
\subsection{ThermosphericDensity}\label{parametrizationAccelerationType:thermosphericDensity}
Estimate the thermospheric density along the orbit using a satllite macro model.
An optional thermospheric model can be used to compute temperature and wind.
The temperature is used to estimate variable drag and lift coefficients, otherwise a constant drag coefficient is used.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "classes/thermosphere/thermosphere.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "classes/miscAccelerations/miscAccelerationsAtmosphericDrag.h"
#include "parametrizationAcceleration.h"

/***** CLASS ***********************************/

/** @brief Estimate the thermospheric density along the orbit.
* @ingroup parametrizationAccelerationGroup
* @see ParametrizationAcceleration */
class ParametrizationAccelerationThermosphericDensity : public ParametrizationAccelerationBase
{
  ThermospherePtr            thermosphere;
  Bool                       useTemperature, useWind;
  Vector3d                   omega;
  ParametrizationTemporalPtr temporal;
  Bool                       perArc;

public:
  ParametrizationAccelerationThermosphericDensity(Config &config);

  Bool isPerArc() const override {return perArc;}
  Bool setInterval(const Time &timeStart, const Time &timeEnd) override {return temporal->setInterval(timeStart, timeEnd, perArc);}
  UInt parameterCount() const override {return temporal->parameterCount();}
  void parameterName(std::vector<ParameterName> &name) const override;
  void compute(SatelliteModelPtr satellite, const Time &time, const Vector3d &position, const Vector3d &velocity,
               const Rotary3d &rotSat, const Rotary3d &rotEarth, EphemeridesPtr ephemerides, MatrixSliceRef A) override;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline ParametrizationAccelerationThermosphericDensity::ParametrizationAccelerationThermosphericDensity(Config &config)
{
  try
  {
    Double angleVelocity;

    readConfig(config, "thermosphere",        thermosphere,   Config::OPTIONAL, "",  "for wind and temperature");
    readConfig(config, "earthRotation",       angleVelocity,  Config::DEFAULT,  "7.29211585531e-5", "[rad/s]");
    readConfig(config, "considerTemperature", useTemperature, Config::DEFAULT,  "1", "compute drag and lift, otherwise simple drag coefficient is used");
    readConfig(config, "considerWind",        useWind,        Config::DEFAULT,  "1", "");
    readConfig(config, "temporalDensity",     temporal,       Config::MUSTSET,  "",  "parameters along orbit");
    readConfig(config, "perArc",              perArc,         Config::DEFAULT,  "0", "");
    if(isCreateSchema(config)) return;

    omega = Vector3d(0, 0, angleVelocity);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void ParametrizationAccelerationThermosphericDensity::parameterName(std::vector<ParameterName> &name) const
{
  try
  {
    temporal->parameterName({ParameterName("satellite", "density")}, name);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void ParametrizationAccelerationThermosphericDensity::compute(SatelliteModelPtr satellite, const Time &time, const Vector3d &position, const Vector3d &velocity,
                                                                     const Rotary3d &rotSat, const Rotary3d &rotEarth, EphemeridesPtr /*ephemerides*/, MatrixSliceRef A)
{
  try
  {
    if(!satellite)
      throw(Exception("No satellite model given"));

    Double   temperature=0;
    Vector3d wind;
    if(thermosphere)
    {
      Double density;
      thermosphere->state(time, rotEarth.rotate(position), density, temperature, wind);
      if(!useTemperature) temperature = 0;
      if(!useWind)        wind = Vector3d();
    }

    // direction and speed of thermosphere relative to satellite in SRF
    Vector3d direction = rotSat.inverseRotate(rotEarth.inverseRotate(wind) + crossProduct(omega, position) - velocity);
    const Double   v   = direction.normalize();
    const Vector3d acc = (1./satellite->mass) * MiscAccelerationsAtmosphericDrag::force(satellite, direction, v, 1/*density*/, temperature);
    temporal->designMatrix(time, rotEarth.rotate(rotSat.rotate(acc)).vector(), A);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}


/***********************************************/

#endif
