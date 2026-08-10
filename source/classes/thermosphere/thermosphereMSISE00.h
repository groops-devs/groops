/***********************************************/
/** @file thermosphereMSISE00.h @brief NRLMSISE-00 thermosphere. */
/***********************************************/
#ifndef __GROOPS_THERMOSPHEREMSISE00__
#define __GROOPS_THERMOSPHEREMSISE00__

#ifdef DOCSTRING_Thermosphere
static const char *docstringThermosphereMSISE00 = R"(
\subsection{MSISE00}
Thermosphere parameters from the NRLMSISE-00 model. The input file uses the
same nine-column layout as \configClass{thermosphere:nrlmsis2}{thermosphereType}:
81-day mean F10.7, daily F10.7, daily Ap, and the six storm-time Ap values.
)";
#endif

#include "base/planets.h"
#include "external/msise00/msise00.h"
#include "classes/thermosphere/thermosphere.h"

class ThermosphereMSISE00 : public Thermosphere
{
  MiscValuesArc msisData;
  Bool interpolateSpaceWeatherIndices;

public:
  inline ThermosphereMSISE00(Config &config);
  inline void state(const Time &time, const Vector3d &position, Double &density,
                    Double &temperature, Vector3d &velocity) const override;
};

inline ThermosphereMSISE00::ThermosphereMSISE00(Config &config)
{
  try
  {
    FileName fileNameMsis, fileNameMagnetic3hAp;
    readConfig(config, "inputfileMsis", fileNameMsis, Config::MUSTSET,
               "{groopsDataDir}/thermosphere/msise00/inputMSIS.txt",
               "F10.7 and Ap input in the NRLMSIS nine-column format");
    readConfig(config, "interpolateSpaceWeatherIndices", interpolateSpaceWeatherIndices,
               Config::OPTIONAL, "0", "linearly interpolate F10.7 and Ap between input epochs");
    readConfig(config, "inputfileMagnetic3hAp", fileNameMagnetic3hAp, Config::OPTIONAL,
               "{groopsDataDir}/thermosphere/hwm14/apActivity.txt", "indices for wind model");
    readConfig(config, "hwm14DataDirectory", fileNameHwm14Path, Config::OPTIONAL,
               "{groopsDataDir}/thermosphere/hwm14",
               "directory containing dwm07b104i.dat, gd2qd.dat, hwm123114.bin");
    if(isCreateSchema(config)) return;

    msisData = InstrumentFile::read(fileNameMsis);
    magnetic3hAp = InstrumentFile::read(fileNameMagnetic3hAp);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

inline void ThermosphereMSISE00::state(const Time &time, const Vector3d &position,
                                        Double &density, Double &temperature,
                                        Vector3d &velocity) const
{
  try
  {
    Ellipsoid ellipsoid;
    Angle lon, lat;
    Double height;
    ellipsoid(position, lon, lat, height);

    const Time timeUt = timeGPS2UTC(time);
    const Vector index = getIndices(msisData, time, interpolateSpaceWeatherIndices);
    const Double dailyF107 = getIndices(msisData, time-mjd2time(1),
                                        interpolateSpaceWeatherIndices)(1);
    const Double averageF107 = index(0);
    const F77Double aps[7] = {index(2), index(3), index(4), index(5),
                              index(6), index(7), index(8)};
    const F77Int doy = static_cast<F77Int>(timeUt.dayOfYear());
    const F77Double sec = timeUt.mjdMod()*86400.;
    const F77Double altitude = height*1e-3;
    const F77Double latitude = lat*RAD2DEG;
    const F77Double longitude = lon*RAD2DEG;
    const F77Double localTime = std::fmod(sec/3600.+longitude/15.+24., 24.);
    F77Double exosphericTemperature;

    msise00CalcWrapper(doy, sec, altitude, latitude, longitude, localTime,
                       averageF107, dailyF107, aps, density, temperature,
                       exosphericTemperature);
    velocity = wind(time, position);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

#endif
