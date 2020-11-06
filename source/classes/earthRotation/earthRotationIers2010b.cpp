/***********************************************/
/**
* @file earthRotationIers2010b.cpp
*
* @brief According to IERS2010 conventions (with HF EOP model).
* @see EarthRotation
*
* @author Torsten Mayer-Guerr
* @date 2019-05-15
*
*/
/***********************************************/

#include "base/import.h"
#include "base/polynomial.h"
#include "external/iers/iers.h"
#include "config/config.h"
#include "files/fileEarthOrientationParameter.h"
#include "files/fileDoodsonEarthOrientationParameter.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/earthRotation/earthRotationIers2010b.h"
#ifndef NOLIB_ERFA
#include <erfa.h>
#endif

/***********************************************/

EarthRotationIers2010b::EarthRotationIers2010b(Config &config)
{
  try
  {
    FileName eopName, doodsonEopName;
    readConfig(config, "inputfileEOP",        eopName,        Config::OPTIONAL,  "{groopsDataDir}/earthRotation/EOP_14C04_IAU2000.txt", "");
    readConfig(config, "inputfileDoodsonEOP", doodsonEopName, Config::OPTIONAL,  "{groopsDataDir}/earthRotation/doodsonEOP_desai_jgrb51665-sup-0002-ds01.txt", "");
    if(isCreateSchema(config)) return;

#ifdef NOLIB_ERFA
    throw(Exception("Compiled without ERFA library"));
#endif

    // read Earth Orientation Parameter (EOP)
    // --------------------------------------
    if(!eopName.empty())
    {
      readFileEarthOrientationParameter(eopName, EOP);
      times.resize(EOP.rows());
      for(UInt i=0; i<times.size(); i++)
        times.at(i) = mjd2time(EOP(i,0));
      EOP = EOP.column(1, 6); // remove mjd

      // UT1-UTC => UT1-GPS (avoid leap seconds jumps for interpolation)
      for(UInt i=0; i<EOP.rows(); i++)
        EOP(i,2) -= (timeUTC2GPS(times.at(i))-times.at(i)).seconds();

      EOP.column(0) *= DEG2RAD/3600; // xp
      EOP.column(1) *= DEG2RAD/3600; // yp
      EOP.column(4) *= DEG2RAD/3600; // dX
      EOP.column(5) *= DEG2RAD/3600; // dY

      polynomial.init(3);
    }

    // read high frequent EOP
    // ----------------------
    if(!doodsonEopName.empty())
    {
      readFileDoodsonEarthOrientationParameter(doodsonEopName, doodsonEop);
      doodsonEop.coeff.column(0, 4) *= DEG2RAD/3600; // xp,yp   [arcsec -> radiant]
      doodsonMatrix = Doodson::matrix(doodsonEop.doodson); // Matrix with Doodson multiplicators
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void EarthRotationIers2010b::earthOrientationParameter(const Time &timeGPS, Double &xp, Double &yp, Double &sp, Double &deltaUT, Double &LOD, Double &X, Double &Y, Double &S) const
{
  try
  {
#ifdef NOLIB_ERFA
    throw(Exception("Compiled without ERFA library"));
#else
    // interpolate EOP file
    // --------------------
    xp = yp = deltaUT = LOD = 0;
    Double dX=0, dY=0;
    if(EOP.size())
    {
      const Time timeUTC = timeGPS2UTC(timeGPS);
      if((timeUTC<times.at(0)) || (timeUTC>times.back()))
        throw(Exception("No EOPs available: "+timeGPS.dateTimeStr()));
      Matrix eop = polynomial.interpolate({timeUTC}, times, EOP, 1);
      xp      = eop(0,0);
      yp      = eop(0,1);
      deltaUT = eop(0,2) + (timeGPS-timeUTC).seconds();
      LOD     = eop(0,3);
      dX      = eop(0,4);
      dY      = eop(0,5);
    }

    // Models
    // ------
    // diurnal and semidiurnal variations in EOP (x,y,UT1) from ocean tides
    if(doodsonEop.doodson.size())
    {
      const Vector thetaf = doodsonMatrix * Doodson::arguments(timeGPS);
      Matrix cosSin(thetaf.rows(),2);
      for(UInt i=0; i<thetaf.rows(); i++)
      {
        cosSin(i, 0) = std::cos(thetaf(i));
        cosSin(i, 1) = std::sin(thetaf(i));
      }
      xp      += inner(cosSin, doodsonEop.coeff.column(0, 2));
      yp      += inner(cosSin, doodsonEop.coeff.column(2, 2));
      deltaUT += inner(cosSin, doodsonEop.coeff.column(4, 2));
      LOD     += inner(cosSin, doodsonEop.coeff.column(6, 2));
    }

    const Double mjdUTC = timeGPS2UTC(timeGPS).mjd();
    Double pm[2];
    pmsdnut2(mjdUTC, pm);
    xp += pm[0]*1e-6*DEG2RAD/3600;
    yp += pm[1]*1e-6*DEG2RAD/3600;

    Double dut1, dlod;
    utlibr(mjdUTC, dut1, dlod);
    deltaUT += dut1*1e-6;
    LOD     += dlod*1e-6;

    const Time timeTT = timeGPS2TT(timeGPS);
    sp = eraSp00(2400000.5+timeTT.mjdInt(), timeTT.mjdMod());

    // precession & nutation
    X = Y = S =0;
    eraXy06   (2400000.5+timeTT.mjdInt(), timeTT.mjdMod(), &X, &Y);
    S = eraS06(2400000.5+timeTT.mjdInt(), timeTT.mjdMod(), X, Y);
    X += dX;
    Y += dY;
#endif
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

