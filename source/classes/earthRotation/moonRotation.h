/***********************************************/
/**
* @file moonRotation.h
*
* @brief Moon rotation using libration angles.
* @see EarthRotation
*
* @author Beate Klinger
* @date 2013-01-01
*
*/
/***********************************************/

#ifndef __GROOPS_MOONROTATIONLIB__
#define __GROOPS_MOONROTATIONLIB__

// Latex documentation
#ifdef DOCSTRING_EarthRotation
static const char *docstringMoonRotation = R"(
\subsection{MoonRotation}
This class realizes the transformation between the moon-fixed system
(Principal Axis System (PA) or Mean Earth System (ME))
and the ICRS according to the JPL ephemeris file.
)";
#endif

/***********************************************/

#include "files/fileEphemerides.h"
#include "classes/earthRotation/earthRotation.h"

/***** CLASS ***********************************/

/** @brief Moon rotation using libration angles.
* @ingroup earthRotationGroup
* @see EarthRotation */
class MoonRotation : public EarthRotation
{
  FileName                  fileName;
  mutable InFileEphemerides file;
  Bool                      useME;

  Rotary3d calcApproxLib(const Time &timeGPS) const;

public:
  MoonRotation(Config &config);

  Rotary3d rotaryMatrix(const Time &timeGPS) const;
};

/***********************************************/
/***********************************************/

MoonRotation::MoonRotation(Config &config)
{
  try
  {
    std::string choice;
    useME = FALSE;

    readConfig(config, "inputfileEphemerides", fileName, Config::OPTIONAL, "{groopsDataDir}/tides/ephemerides_JPL_DE432.dat", "librations");
    if(readConfigChoice(config, "moonfixedSystem", choice, Config::MUSTSET, "", ""))
    {
      if(readConfigChoiceElement(config, "PA", choice, "Principal Axis System")) {useME = FALSE;}
      if(readConfigChoiceElement(config, "ME", choice, "Mean Earth System"))     {useME = TRUE;}
      endChoice(config);
    }
    if(isCreateSchema(config)) return;

    file.open(fileName);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Rotary3d MoonRotation::rotaryMatrix(const Time &timeGPS) const
{
  try
  {
    Rotary3d rot;
    if(fileName.empty())
    {
      rot = calcApproxLib(timeGPS);
    }
    else
    {
      const Matrix A = file.interpolate(timeGPS, InFileEphemerides::LIBRATION);
      rot = rotaryZ(Angle(A(2,0))) * rotaryX(Angle(A(1,0))) * rotaryZ(Angle(A(0,0)));
    }

    if(useME)
      rot = rotaryX(Angle((-0.3/3600)*DEG2RAD)) * rotaryY(Angle((-78.56/3600)*DEG2RAD)) * rotaryZ(Angle((-67.92/3600)*DEG2RAD)) * rot;
    return rot;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Rotary3d MoonRotation::calcApproxLib (const Time &timeGPS) const
{
  const Double d = timeGPS2TT(timeGPS).mjd()-J2000; // interval in days from standard epoch (J2000.0)
  const Double T = d/36525;             // interval in Julian centuries (36525 days) from standard epoch

  const Double E1  = (125.045*DEG2RAD) - ( 0.0529921*DEG2RAD)*d;
  const Double E2  = (250.089*DEG2RAD) - ( 0.1059842*DEG2RAD)*d;
  const Double E3  = (260.008*DEG2RAD) + (13.0120009*DEG2RAD)*d;
  const Double E4  = (176.625*DEG2RAD) + (13.3407154*DEG2RAD)*d;
  const Double E5  = (357.529*DEG2RAD) + ( 0.9856003*DEG2RAD)*d;
  const Double E6  = (311.589*DEG2RAD) + (26.4057084*DEG2RAD)*d;
  const Double E7  = (134.963*DEG2RAD) + (13.0649930*DEG2RAD)*d;
  const Double E8  = (276.617*DEG2RAD) + ( 0.3287146*DEG2RAD)*d;
  const Double E9  = ( 34.226*DEG2RAD) + ( 1.7484877*DEG2RAD)*d;
  const Double E10 = ( 15.134*DEG2RAD) - ( 0.1589763*DEG2RAD)*d;
  const Double E11 = (119.743*DEG2RAD) + ( 0.0036096*DEG2RAD)*d;
  const Double E12 = (239.961*DEG2RAD) + ( 0.1643573*DEG2RAD)*d;
  const Double E13 = ( 25.053*DEG2RAD) + (12.9590088*DEG2RAD)*d;

  // ICRF equatoral coordinates at epoch J2000.0
  const Double alpha0 = (269.9949*DEG2RAD) + (0.0031*DEG2RAD)*T - (3.8787*DEG2RAD)*sin(E1)
                      - (0.1204*DEG2RAD)*sin(E2) + (0.0700*DEG2RAD)*sin(E3) - (0.0172*DEG2RAD)*sin(E4) + (0.0072*DEG2RAD)*sin(E6)
                      - (0.0052*DEG2RAD)*sin(E10) + (0.0043*DEG2RAD)*sin(E13);

  const Double delta0 = (66.5392*DEG2RAD) + (0.0130*DEG2RAD)*T + (1.5419*DEG2RAD)*cos(E1) + (0.0239*DEG2RAD)*cos(E2) - (0.0278*DEG2RAD)*cos(E3)
                      + (0.0068*DEG2RAD)*cos(E4) - (0.0029*DEG2RAD)*cos(E6) + (0.0009*DEG2RAD)*cos(E7) + (0.0008*DEG2RAD)*cos(E10) - (0.0009*DEG2RAD)*cos(E13);


  const Double W = (38.3213*DEG2RAD) + (13.17635815*DEG2RAD)*d - (1.4*pow(10,-12)*DEG2RAD)*pow(d,2) + (3.5610*DEG2RAD)*sin(E1) + (0.1208*DEG2RAD)*sin(E2)
                 - (0.0642*DEG2RAD)*sin(E3) + (0.0158*DEG2RAD)*sin(E4) + (0.0252*DEG2RAD)*sin(E5) - (0.0066*DEG2RAD)*sin(E6) - (0.0047*DEG2RAD)*sin(E7)
                 - (0.0046*DEG2RAD)*sin(E8) + (0.0028*DEG2RAD)*sin(E9) + (0.0052*DEG2RAD)*sin(E10) + (0.0040*DEG2RAD)*sin(E11) + (0.0019*DEG2RAD)*sin(E12) - (0.0044*DEG2RAD)*sin(E13);


  return rotaryZ(Angle(W)) * rotaryX(Angle(PI/2-delta0)) * rotaryZ(Angle(alpha0+PI/2));  // Rotationsmatrix: ICRF >> ME
}

/***********************************************/

#endif /* __GROOPS_GRAILMOONROTATIONLIB__ */
