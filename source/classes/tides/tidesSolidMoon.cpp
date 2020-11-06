/***********************************************/
/**
* @file tidesSolidMoon.cpp
*
* @brief Solid moon tides.
* @see Tides
*
* @author Beate Klinger
* @date 2013-xx-xx
*
*/
/***********************************************/

#include "base/import.h"
#include "base/sphericalHarmonics.h"
#include "config/config.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/tides/tidesSolidMoon.h"

/***********************************************/

TidesSolidMoon::TidesSolidMoon(Config &config)
{
  try
  {
    readConfig(config, "k20",    kReal20, Config::DEFAULT,  "0.0213", "");
    readConfig(config, "k30",    kReal30, Config::DEFAULT,  "0.0",    "");
    readConfig(config, "factor", factor,  Config::DEFAULT,  "1.0",    "the result is multplied by this factor, set -1 to substract the field");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/************************************************/

// Geopotential (IERS-Convetions 2010) Step 1
// ------------------------------------------------
void TidesSolidMoon::moonCoefficients1(Double GM_third, const Vector3d &third, Matrix &cnm, Matrix &snm) const
{
  // Kugelflaechenfunktionen berechnen
  Matrix Cnm, Snm;
  SphericalHarmonics::CnmSnm(1./R_Moon * third, 3, Cnm, Snm);

  Double factor1 = GM_third/GM_Moon;
  Double kReal = 0.0;

  // Formel (1)
  for(UInt n=2; n<=3; n++)
    for(UInt m=0; m<=0; m++)
    {
      if (n == 2)
        kReal = kReal20;
      if (n == 3)
        kReal = kReal30;

      // nur Realteil berücksichtigt, kein Imaginärteil
      cnm(n,m) += factor1/(2.*n+1.) * (kReal * Cnm(n,m));
      snm(n,m) += factor1/(2.*n+1.) * (kReal * Snm(n,m));
    }
}

/***********************************************/

SphericalHarmonics TidesSolidMoon::sphericalHarmonics(const Time &time, const Rotary3d &rotMoon, EarthRotationPtr /*rotation*/, EphemeridesPtr ephemerides, UInt maxDegree, UInt minDegree, Double GM, Double R) const
{
  try
  {
    Matrix cnm(5, Matrix::TRIANGULAR, Matrix::LOWER);
    Matrix snm(5, Matrix::TRIANGULAR, Matrix::LOWER);

    const Vector3d posMoon = ephemerides->position(time, Ephemerides::MOON);
    moonCoefficients1(GM_Sun,   rotMoon.rotate(ephemerides->position(time, Ephemerides::SUN)   - posMoon), cnm, snm);
    moonCoefficients1(GM_Earth, rotMoon.rotate(ephemerides->position(time, Ephemerides::EARTH) - posMoon), cnm, snm);

    return SphericalHarmonics(GM_Moon, R_Moon, factor*cnm, factor*snm).get(maxDegree, minDegree, GM, R);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
