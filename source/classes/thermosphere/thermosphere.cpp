/***********************************************/
/**
* @file thermosphere.cpp
*
* @brief Density, temperature and velocity.
*
* @author Torsten Mayer-Guerr
* @author Andreas Kvas
* @date 2020-02-20
*
*/
/***********************************************/

#define DOCSTRING_Thermosphere

#include "base/import.h"
#include "external/hwm/hwm.h"
#include "config/configRegister.h"
#include "classes/thermosphere/thermosphereJB2008.h"
#include "classes/thermosphere/thermosphere.h"


/***********************************************/

GROOPS_REGISTER_CLASS(Thermosphere, "thermosphereType",
                      ThermosphereJB2008)

GROOPS_READCONFIG_CLASS(Thermosphere, "thermosphereType")

/***********************************************/

ThermospherePtr Thermosphere::create(Config &config, const std::string &name)
{
  try
  {
    ThermospherePtr thermosphere;
    std::string choice;

    readConfigChoice(config, name, choice, Config::MUSTSET, "", "density, temperature and velocity");
    if(readConfigChoiceElement(config, "jb2008",  choice, "Jacchia-Bowman 2008 Empirical Thermospheric Density Model"))
      thermosphere = ThermospherePtr(new ThermosphereJB2008(config));
    endChoice(config);

    return thermosphere;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector Thermosphere::getIndices(const MiscValuesArc &arc, const Time &time, Bool interpolate)
{
  try
  {
    const Time timeUt = timeGPS2UTC(time);
    auto iter = std::upper_bound(arc.begin(), arc.end(), timeUt, [](const Time &time, auto &epoch) {return time < epoch.time;});
    if((iter == arc.begin()) || (iter == arc.end()))
      throw(Exception("Cannot find index for "+time.dateTimeStr()));
    iter--;
    if(!interpolate)
      return iter->values;
    const Double t = (timeUt - iter->time).seconds()/((iter+1)->time - iter->time).seconds();
    return (1-t) * iter->values + t * (iter+1)->values;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector3d Thermosphere::wind(const Time &time, const Vector3d &position) const
{
  try
  {
    if(fileNameHwm14Path.empty())
      return Vector3d();

    Double ap2 = -1;
    if(magnetic3hAp.size())
      ap2 = getIndices(magnetic3hAp, time, FALSE)(0);

    Ellipsoid ellipsoid;
    Angle     lon, lat;
    Double    height;
    ellipsoid(position, lon, lat, height);

    UInt   year, month, day, hour, min;
    Double seconds;
    const Time timeUt = timeGPS2UTC(time);
    timeUt.date(year, month, day, hour, min, seconds);
    F77Int yyddd = ((year<2000) ? (year-1900) : (year-2000))*1000 + timeUt.dayOfYear();

#ifdef _WIN32
    _putenv_s("HWMPATH", fileNameHwm14Path.c_str());
#else
    setenv("HWMPATH", fileNameHwm14Path.c_str(), 1);
#endif
    F77Float outf[2];
    const F77Float ap[] = {0.0, static_cast<F77Float>(ap2)};
    hwm(yyddd, 3600*hour+60*min+seconds, height*0.001, lat*RAD2DEG, std::fmod(lon*RAD2DEG+360, 360), 0, 0., 0., ap, outf);
#ifndef _WIN32
    unsetenv("HWMPATH");
#endif

    return localNorthEastUp(position, ellipsoid).transform(Vector3d(outf[0], outf[1], 0));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
