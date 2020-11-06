/***********************************************/
/**
* @file ephemeridesJpl.h
*
* @brief Ephemerides from JPL.
*
* @author Torsten Mayer-Guerr
* @date 2020-06-07
*
*/
/***********************************************/

#ifndef __GROOPS_EPHEMERIDESJPL__
#define __GROOPS_EPHEMERIDESJPL__

// Latex documentation
#ifdef DOCSTRING_Ephemerides
static const char *docstringEphemeridesJpl = R"(
\section{JPL}\label{ephemeridesType:jpl}
Using \verb|DExxx| ephemerides from NASA Jet Propulsion Laboratory (JPL).
)";
#endif

/***********************************************/

#include "files/fileEphemerides.h"
#include "classes/ephemerides/ephemerides.h"

/***** CLASS ***********************************/

/** @brief Ephemerides from JPL.
* @ingroup ephemeridesGroup
* @see Ephemerides */
class EphemeridesJpl : public Ephemerides
{
  InFileEphemerides file;
  Planet origin_;

public:
  EphemeridesJpl(Config &config);

  Vector3d position(const Time &timeGPS, Planet planet) override;
  void     ephemeris(const Time &timeGPS, Planet planet, Vector3d &position, Vector3d &velocity) override;
  Planet   origin() const override {return origin_;}
};

/***********************************************/

inline EphemeridesJpl::EphemeridesJpl(Config &config)
{
  try
  {
    FileName fileName;

    readConfig(config, "inputfileEphemerides", fileName, Config::MUSTSET, "{groopsDataDir}/tides/ephemerides_JPL_DE432.dat", "");
    readConfig(config, "origin",               origin_,  Config::DEFAULT, "earth", "center of coordinate system");
    if(isCreateSchema(config)) return;

    file.open(fileName);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Vector3d EphemeridesJpl::position(const Time &timeGPS, Planet planet)
{
  try
  {
    Vector3d position, velocity;
    file.ephemeris(timeGPS, static_cast<InFileEphemerides::Planet>(planet), static_cast<InFileEphemerides::Planet>(origin_), position, velocity);
    return position;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void EphemeridesJpl::ephemeris(const Time &timeGPS, Planet planet, Vector3d &position, Vector3d &velocity)
{
  try
  {
    file.ephemeris(timeGPS, static_cast<InFileEphemerides::Planet>(planet), static_cast<InFileEphemerides::Planet>(origin_), position, velocity);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/************************************************/

#endif
