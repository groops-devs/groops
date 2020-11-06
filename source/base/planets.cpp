/***********************************************/
/**
* @file planets.cpp
*
* @brief Planet and Earth fundamentals.
*
* @author Torsten Mayer-Guerr
* @date 2005-04-26
*
*/
/***********************************************/

#include "base/importStd.h"
#include "base/matrix.h"
#include "base/time.h"
#include "base/vector3d.h"
#include "base/transform3d.h"
#include "base/ellipsoid.h"
#include "base/planets.h"

/************************************************/

// Reference: The Astronomical Almanac, page C24.
Vector3d Planets::positionSun(const Time &timeGPS)
{
  const Double t = (timeGPS2TT(timeGPS)-mjd2time(J2000)).mjd();
  // mean anomaly
  const Double g = DEG2RAD * (357.528 + 0.9856003*t);
  // solar ecliptic longitude
  const Double lambda = DEG2RAD * (280.464 + 0.9856090*t + 1.915*std::sin(g) + 0.02*std::sin(2*g));
  // obliquity of the ecliptic
  const Double eps = DEG2RAD * (23.439 - 0.0000004*t);
  // distance
  const Double r = 149597870700 * (1.000142 - 0.01671*std::cos(g) - 0.00014*std::cos(2*g));
  // position
  return Vector3d(r*std::cos(lambda), r*std::sin(lambda)*std::cos(eps), r*std::sin(lambda)*std::sin(eps));
}

/************************************************/

// Fundamental Arguments of nutation theory
// from IERS Conventions 2003
Vector Planets::fundamentals(const Time &timeGPS)
{
  const Double T = timeGPS2JC(timeGPS);
  Vector F(5);
  F(0) = DEG2RAD/3600*( 485868.249036+T*(1717915923.2178+T*( 31.8792+T*( 0.051635+T*(-0.00024470)))));
  F(1) = DEG2RAD/3600*(1287104.793048+T*( 129596581.0481+T*( -0.5532+T*( 0.000136+T*(-0.00001149)))));
  F(2) = DEG2RAD/3600*( 335779.526232+T*(1739527262.8478+T*(-12.7512+T*(-0.001037+T*( 0.00000417)))));
  F(3) = DEG2RAD/3600*(1072260.703692+T*(1602961601.2090+T*( -6.3706+T*( 0.006593+T*(-0.00003169)))));
  F(4) = DEG2RAD/3600*( 450160.398036+T*( - 6962890.5431+T*(  7.4722+T*( 0.007702+T*(-0.00005939)))));
  return F;
}

/************************************************/

Double Planets::gmst(const Time &timeUT1)
{
  Double Tu0 = (timeUT1.mjdInt()-51544.5)/36525.0;

  Double GMST0 = (6.0/24 + 41.0/(24*60) + 50.54841/(24*60*60))
               + (8640184.812866/(24*60*60))*Tu0
               + (0.093104/(24*60*60))*Tu0*Tu0
               + (-6.2e-6/(24*60*60))*Tu0*Tu0*Tu0;
  Double r     = 1.002737909350795 + 5.9006e-11*Tu0 - 5.9e-15*Tu0*Tu0;
  return fmod(2*PI*(GMST0 + r * timeUT1.mjdMod()), 2*PI);
}

/************************************************/

Double Planets::ERA(const Time &timeUT1)
{
  const Time T = timeUT1-mjd2time(J2000);
  return fmod(2*PI*(0.7790572732640 + T.mjdMod() + 0.00273781191135448*T.mjd()), 2*PI);
}

/************************************************/

Rotary3d Planets::celestial2TerrestrialFrame(const Time &timeGPS)
{
  return rotaryZ(Angle(Planets::gmst(timeGPS2UTC(timeGPS))));
}

/***********************************************/

Double Planets::normalGravity(const Vector3d &p)
{
  const Double ga = 9.7803267715;
  const Double gb = 9.8321863685;
  const Double m  = 0.00344978600308;
  const Double f  = DEFAULT_GRS80_f;
  const Double a  = DEFAULT_GRS80_a;
  const Double b  = a*(1-1/f);

  Ellipsoid ellipsoid(DEFAULT_GRS80_a, DEFAULT_GRS80_f);
  Angle  L,B;
  Double h;
  ellipsoid(p, L,B,h);

  const Double cos2 = std::pow(std::cos(B),2);
  const Double sin2 = std::pow(std::sin(B),2);
  const Double gamma0 = (a*ga*cos2+b*gb*sin2)/std::sqrt(a*a*cos2+b*b*sin2);
  return gamma0 -2*ga/a*(1+1/f+m+(-3/f+5*m/2)*sin2)*h+3*ga/a/a*h*h;
}

/************************************************/
