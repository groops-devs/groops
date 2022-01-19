/***********************************************/
/**
* @file thermosphereJB2008.h
*
* @brief Density, temperature and velocity.
*
* @author Torsten Mayer-Guerr
* @date 2020-02-20
*
*/
/***********************************************/

#ifndef __GROOPS_THERMOSPHEREJB2008__
#define __GROOPS_THERMOSPHEREJB2008__

// Latex documentation
#ifdef DOCSTRING_Thermosphere
static const char *docstringThermosphereJB2008 = R"(
\subsection{JB2008}
Thermosphere parameters from the JB2008 model:

Bowman, B. R., Tobiska, W. K., Marcos, F. A., Huang, C. Y., Lin, C. S., Burke, W. J. (2008).
A new empirical thermospheric density model JB2008 using new solar and geomagnetic indices.
 In AIAA/AAS Astrodynamics Specialist Conference and Exhibit. \url{https://doi.org/10.2514/6.2008-6438}
)";
#endif

/***********************************************/

#include "base/planets.h"
#include "external/jb2008/jb2008.h"
#include "inputOutput/file.h"
#include "classes/thermosphere/thermosphere.h"

/***** CLASS ***********************************/

/** @brief Density, temperature and velocity.
* @ingroup thermosphereGroup
* @see Thermosphere */
class ThermosphereJB2008 : public Thermosphere
{
  MiscValuesArc solarFSMY;
  MiscValuesArc dtc;

public:
  inline ThermosphereJB2008(Config &config);

  inline void state(const Time &time, const Vector3d &position, Double &density, Double &temperature, Vector3d &velocity) const override;
};

/***********************************************/

inline ThermosphereJB2008::ThermosphereJB2008(Config &config)
{
  try
  {
    FileName fileNameSolfsmy, fileNameDtc, fileNameMagnetic3hAp;

    readConfig(config, "inputfileSolfsmy",      fileNameSolfsmy,      Config::MUSTSET,  "{groopsDataDir}/thermosphere/jb2008/SOLFSMY.TXT",   "solar indices");
    readConfig(config, "inputfileDtc",          fileNameDtc,          Config::MUSTSET,  "{groopsDataDir}/thermosphere/jb2008/DTCFILE.TXT",   "");
    readConfig(config, "inputfileMagnetic3hAp", fileNameMagnetic3hAp, Config::OPTIONAL, "{groopsDataDir}/thermosphere/hwm14/apActivity.txt", "indicies for wind model");
    readConfig(config, "hwm14DataDirectory",    fileNameHwm14Path,    Config::OPTIONAL, "{groopsDataDir}/thermosphere/hwm14",                "directory containing dwm07b104i.dat, gd2qd.dat, hwm123114.bin");
    if(isCreateSchema(config)) return;

#ifdef GROOPS_DISABLE_JB2008
    throw(Exception("Compiled without JB2008 sources"));
#endif

    if(!fileNameSolfsmy.empty())
    {
      InFile file(fileNameSolfsmy);
      std::string line;
      for(UInt i=0; i<4; i++)
        std::getline(file, line);
      while(std::getline(file, line))
      {
        std::stringstream ss(line);
        UInt        year, day;
        Double      jd, f, fb, s, sb, m, mb, y, yb;
        std::string flags;
        ss>>year>>day>>jd>>f>>fb>>s>>sb>>m>>mb>>y>>yb>>flags;

        MiscValuesEpoch epoch(8);
        epoch.time = mjd2time(jd-2400001);
        epoch.values = {f, fb, s, sb, m, mb, y, yb};
        solarFSMY.push_back(epoch);
      }
    }

    if(!fileNameDtc.empty())
    {
      InFile file(fileNameDtc);
      std::string line;
      while(std::getline(file, line))
      {
        std::stringstream ss(line);
        std::string tag;
        UInt        year, day;
        ss>>tag>>year>>day;
        for(UInt i=0; i<24; i++)
        {
          MiscValuesEpoch epoch(1);
          epoch.time = date2time(year, 1, 1, i, 0, 0.) + mjd2time(day-1.);
          ss>>epoch.values(0);
          dtc.push_back(epoch);
        }
      }
    }

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

inline void ThermosphereJB2008::state(const Time &time, const Vector3d &position, Double &density, Double &temperature, Vector3d &velocity) const
{
  try
  {
#ifndef GROOPS_DISABLE_JB2008
    Ellipsoid ellipsoid;
    Angle     lon, lat;
    Double    height;
    ellipsoid(position, lon, lat, height);

    // use 5 day lag for y10 for jb2008
    const Vector index5 = getIndices(solarFSMY, time - mjd2time(5), FALSE);
    const Double Y10  = index5(6);
    const Double Y10B = index5(7);

    // use 2 day lag for m10 for jb2008
    const Vector index2 = getIndices(solarFSMY, time - mjd2time(2), FALSE);
    const Double M10  = index2(4);
    const Double M10B = index2(5);

    // use 1 day lag for f10 and s10 for jb2008
    const Vector index1 = getIndices(solarFSMY, time - mjd2time(1), FALSE);
    const Double F10  = index1(0);
    const Double F10B = index1(1);
    const Double S10  = index1(2);
    const Double S10B = index1(3);

    // read geomagnetic storm dtc value
    const Double dstdtc = getIndices(dtc, time, TRUE)(0);

    const Vector3d sunPos = Planets::positionSun(time);

    F77Double sun[2] = {sunPos.lambda(), sunPos.phi()};
    F77Double pos[3] = {std::fmod(Double(lon+Planets::gmst(timeGPS2UTC(time)))+2*PI, 2*PI), lat, height*1e-3};
    F77Double temp[2], rho;

    jb2008(time.mjd(), sun, pos, F10, F10B, S10, S10B, M10, M10B, Y10, Y10B, dstdtc, temp, rho);

    density     = rho;
    temperature = temp[1];
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
