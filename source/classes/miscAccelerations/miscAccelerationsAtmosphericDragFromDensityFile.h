/***********************************************/
/**
* @file miscAccelerationsAtmosphericDragFromDensityFile.h
*
* @brief Atmospheric drag from denity along orbit.
* @see MiscAccelerations
*
* @author Torsten Mayer-Guerr
* @date 2022-08-13
*
*/
/***********************************************/

#ifndef __GROOPS_MISCACCELERATIONSATMOSPHERICDRAGFROMDENSITYFILE__
#define __GROOPS_MISCACCELERATIONSATMOSPHERICDRAGFROMDENSITYFILE__

// Latex documentation
#ifdef DOCSTRING_MiscAccelerations
static const char *docstringMiscAccelerationsAtmosphericDragFromDensityFile = R"(
\subsection{AtmosphericDrag}\label{miscAccelerationsType:atmosphericDragFromDensityFile}
Atmospheric drag computed from thermospheric density along the orbit
(\configFile{inputfileDensity}{instrument}, MISCVALUE). The \configClass{thermosphere}{thermosphereType}
is used to to compute temperature and wind.
For further details see \configClass{atmosphericDrag}{miscAccelerationsType:atmosphericDrag}.
)";
#endif

/***********************************************/

#include "files/fileInstrument.h"
#include "classes/thermosphere/thermosphere.h"
#include "classes/miscAccelerations/miscAccelerations.h"

/***** CLASS ***********************************/

/** @brief Atmospheric drag.
 * @ingroup miscAccelerationsGroup
 * @see MiscAccelerations */
class MiscAccelerationsAtmosphericDragFromDensityFile : public MiscAccelerationsBase
{
  MiscValueArc    density;
  ThermospherePtr thermosphere;
  Bool            useTemperature, useWind;
  Vector3d        omega;
  Double          factor;
  UInt            idx;

public:
  MiscAccelerationsAtmosphericDragFromDensityFile(Config &config);

  Vector3d acceleration(SatelliteModelPtr satellite, const Time &time, const Vector3d &position, const Vector3d &velocity,
                        const Rotary3d &rotSat, const Rotary3d &rotEarth, EphemeridesPtr ephemerides) override;
};

/***********************************************/

inline MiscAccelerationsAtmosphericDragFromDensityFile::MiscAccelerationsAtmosphericDragFromDensityFile(Config &config)
{
  try
  {
    FileName fileNameDensity;
    Double   angleVelocity;

    readConfig(config, "inputfileDensity",    fileNameDensity, Config::MUSTSET, "",    "density along orbit, MISCVALUE (kg/m^3)");
    readConfig(config, "thermosphere",        thermosphere,    Config::MUSTSET, "",    "used to compute temperature and wind");
    readConfig(config, "earthRotation",       angleVelocity,   Config::DEFAULT, "7.29211585531e-5", "[rad/s]");
    readConfig(config, "considerTemperature", useTemperature,  Config::DEFAULT, "1",   "compute drag and lift, otherwise simple drag coefficient is used");
    readConfig(config, "considerWind",        useWind,         Config::DEFAULT, "1",   "");
    readConfig(config, "factor",              factor,          Config::DEFAULT, "1.0", "the result is multiplied by this factor");
    if(isCreateSchema(config)) return;

    density = InstrumentFile::read(fileNameDensity);
    omega = Vector3d(0, 0, angleVelocity);
    idx = 0;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Vector3d MiscAccelerationsAtmosphericDragFromDensityFile::acceleration(SatelliteModelPtr satellite, const Time &time,
                                                                              const Vector3d &position, const Vector3d &velocity,
                                                                              const Rotary3d &rotSat, const Rotary3d &rotEarth, EphemeridesPtr /*ephemerides*/)
{
  try
  {
    if(!satellite)
      throw(Exception("No satellite model given"));

    if((time < density.front().time)|| (time > density.back().time))
      throw(Exception("time not given in density file: "+time.dateTimeStr()));

    // find index (interpolations interval)
    if((idx >= density.size()) || (time < density.at(idx).time))
      idx = 0;
    while(time > density.at(idx).time)
      idx++;
    if(time != density.at(idx).time)
      throw(Exception("time not given in density file: "+time.dateTimeStr()));

    Double   densityModel, temperature;
    Vector3d wind;
    thermosphere->state(time, rotEarth.rotate(position), densityModel, temperature, wind);
    if(!useTemperature)
      temperature = 0;
    if(!useWind)
      wind = Vector3d();

    const Vector3d velocityRelativeToThermosphere = velocity - crossProduct(omega, position) - rotEarth.inverseRotate(wind);
    const Vector3d acc = satellite->accelerationDrag(rotSat.inverseRotate(velocityRelativeToThermosphere), density.at(idx).value, temperature);
    return factor * rotEarth.rotate(rotSat.rotate(acc));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
