/***********************************************/
/**
* @file thermosphereNRLMSIS2.h
*
* @brief Density, temperature and velocity.
*
* @author Sandro Krauss
* @date 2021-03-22
*
*/
/***********************************************/

#ifndef __GROOPS_THERMOSPHERENRLMSIS2__
#define __GROOPS_THERMOSPHERENRLMSIS2__

// Latex documentation
#ifdef DOCSTRING_Thermosphere
static const char *docstringThermosphereNRLMSIS2 = R"(
\subsection{NRLMSIS2}
Thermosphere parameters from the NRLMSIS2 model:

Emmert J.D, D.P.Drob, J.M. Picone, et al. (2020), NRLMSIS 2.0: A whole-atmosphere empirical
model of temperature and neutral species densities. Earth and Space Science, Volume 8, 3
\url{https://doi.org/10.1029/2020EA001321}
)";
#endif

/***********************************************/

#include "base/planets.h"
#include "external/nrlmsis2/nrlmsis2.h"
#include "inputOutput/file.h"
#include "classes/thermosphere/thermosphere.h"

/***** CLASS ***********************************/

/** @brief Density, temperature and velocity.
* @ingroup thermosphereGroup
* @see Thermosphere */
class ThermosphereNRLMSIS2 : public Thermosphere
{
  MiscValuesArc msisData;

public:
  inline ThermosphereNRLMSIS2(Config &config);

  inline void state(const Time &time, const Vector3d &position, Double &density, Double &temperature, Vector3d &velocity) const override;
};

/***********************************************/

inline ThermosphereNRLMSIS2::ThermosphereNRLMSIS2(Config &config)
{
  try
  {
    FileName fileNameMsis, fileNameParm, fileNameMagnetic3hAp;

    readConfig(config, "inputfileMsis",            fileNameMsis,         Config::MUSTSET,  "{groopsDataDir}/thermosphere/nrlmsis2/inputMSIS.txt", "input NRLMSIS 2.0");
    readConfig(config, "inputfileModelParameters", fileNameParm,         Config::MUSTSET,  "{groopsDataDir}/thermosphere/nrlmsis2/msis20.parm",   "path to msis20.parm file");
    readConfig(config, "inputfileMagnetic3hAp",    fileNameMagnetic3hAp, Config::OPTIONAL, "{groopsDataDir}/thermosphere/hwm14/apActivity.txt",   "indicies for wind model");
    readConfig(config, "hwm14DataDirectory",       fileNameHwm14Path,    Config::OPTIONAL, "{groopsDataDir}/thermosphere/hwm14",                  "directory containing dwm07b104i.dat, gd2qd.dat, hwm123114.bin");
    if(isCreateSchema(config)) return;

#ifdef GROOPS_DISABLE_NRLMSIS
    throw(Exception("Compiled without NRLMSIS sources"));
#endif

    msisData = InstrumentFile::read(fileNameMsis);
    msisinitWrapper(fileNameParm.directory(), fileNameParm.stripDirectory());

    magnetic3hAp = InstrumentFile::read(fileNameMagnetic3hAp);
#ifdef GROOPS_DISABLE_HWM14
    if(!fileNameHwm14Path.empty())
      logWarningOnce<<"Compiled without HWM14 wind model sources -> thermospheric wind is not calculated"<<Log::endl;
#endif
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void ThermosphereNRLMSIS2::state(const Time &time, const Vector3d &position, Double &density, Double &temperature, Vector3d &velocity) const
{
  try
  {
#ifndef GROOPS_DISABLE_NRLMSIS
    Ellipsoid ellipsoid;
    Angle     lon, lat;
    Double    height;
    ellipsoid(position, lon, lat, height);
    if(lon < 0)
      lon += Angle(2*PI);

    // get data
    const Time     timeUt      =  timeGPS2UTC(time);
    const Vector   index       =  getIndices(msisData, time, FALSE);
    const F77Float dailyF107   =  static_cast<F77Float>(getIndices(msisData, time-mjd2time(1), FALSE)(1));
    const F77Float averageF107 =  static_cast<F77Float>(index(0));
    const F77Float aps[7]      = {static_cast<F77Float>(index(2)),
                                  static_cast<F77Float>(index(3)),
                                  static_cast<F77Float>(index(4)),
                                  static_cast<F77Float>(index(5)),
                                  static_cast<F77Float>(index(6)),
                                  static_cast<F77Float>(index(7)),
                                  static_cast<F77Float>(index(8))};

    // CALL NRLMSIS 2.0 model
    F77Float tempAlt = 0, tempExo = 0, rho[10] = {0.};
    msiscalcWrapper(timeUt.dayOfYear()+timeUt.mjdMod(), time.mjdMod()*86400, height*1e-3, lat*RAD2DEG, lon*RAD2DEG, averageF107, dailyF107, aps, tempAlt, rho, tempExo);
    density     = rho[0];
    temperature = tempAlt;
    velocity    = wind(time, position);
#endif
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
