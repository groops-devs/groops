/***********************************************/
/**
* @file earthRotationIers1996.cpp
*
* @brief According to IERS1996 conventions.
* @see EarthRotation
*
* @author Torsten Mayer-Guerr
* @date 2003-09-20
*
*/
/***********************************************/

#include "base/import.h"
#include "base/polynomial.h"
#include "base/planets.h"
#include "external/iers/iers.h"
#include "config/config.h"
#include "inputOutput/fileArchive.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/earthRotation/earthRotationIers1996.h"

/***********************************************/

EarthRotationIers1996::EarthRotationIers1996(Config &config)
{
  try
  {
    FileName eopName, nutationName;
    readConfig(config, "inputfileEOP",      eopName,      Config::OPTIONAL, "", "");
    readConfig(config, "inputfileNutation", nutationName, Config::MUSTSET,  "{groopsDataDir}/earthRotation/nutationIAU1980.xml", "");
    if(isCreateSchema(config)) return;

#ifdef GROOPS_DISABLE_IERS
    logWarningOnce<<"Compiled without IERS sources -> ocean tidal effects in EOP are not calculated"<<Log::endl;
#endif

    // read Earth Orientation Parameter (EOP)
    // --------------------------------------
    if(!eopName.empty())
    {
      InFileArchive eopFile(eopName, "", FILE_BASE_VERSION);
      UInt eopCount;
      eopFile>>nameValue("count", eopCount);
      std::vector<Time> times(eopCount);
      EOP = Matrix(eopCount, 7);
      Double mjd;
      for(UInt i=0; i<eopCount; i++)
      {
        eopFile>>beginGroup("eop");
        eopFile>>nameValue("mjd", mjd);
        times.at(i) = mjd2time(mjd);
        eopFile>>nameValue("xp",       EOP(i,0));
        eopFile>>nameValue("yp",       EOP(i,1));
        eopFile>>nameValue("deltaUT",  EOP(i,2));
        eopFile>>nameValue("lod",      EOP(i,3));
        eopFile>>nameValue("deltaPsi", EOP(i,4));
        eopFile>>nameValue("deltaEps", EOP(i,5));
        eopFile>>endGroup("eop");
      }

      // UT1-UTC => UT1-GPS (avoid leap seconds jumps for interpolation)
      for(UInt i=0; i<eopCount; i++)
        EOP(i,2) -= (timeUTC2GPS(times.at(i))-times.at(i)).seconds();

      EOP.column(0) *= DEG2RAD/3600; // xp
      EOP.column(1) *= DEG2RAD/3600; // yp
      EOP.column(4) *= DEG2RAD/3600; // dPsi
      EOP.column(5) *= DEG2RAD/3600; // dEps

      polynomial.init(times, 3);
    }

    // Nutationsserie einlesen
    // -----------------------
    InFileArchive nutationFile(nutationName, "", FILE_BASE_VERSION);
    UInt nutationCount;
    nutationFile>>nameValue("count", nutationCount);
    argument  = Matrix(nutationCount, 5);
    psiFactor = Matrix(nutationCount, 2);
    epsFactor = Matrix(nutationCount, 2);

    for(UInt i=0; i<nutationCount; i++)
    {
      nutationFile>>beginGroup("nutation");
      nutationFile>>nameValue("arg1",argument(i,0))>>nameValue("arg2",argument(i,1))>>nameValue("arg3",argument(i,2))>>nameValue("arg4",argument(i,3))>>nameValue("arg5",argument(i,4));
      nutationFile>>nameValue("A1",psiFactor(i,0))>>nameValue("A2",psiFactor(i,1));
      nutationFile>>nameValue("B1",epsFactor(i,0))>>nameValue("B2",epsFactor(i,1));
      nutationFile>>endGroup("nutation");
    }

    psiFactor *= 1e-4*DEG2RAD/3600;
    epsFactor *= 1e-4*DEG2RAD/3600;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Rotary3d EarthRotationIers1996::rotaryMatrix(const Time &timeGPS) const
{
  // Erdrotationsparameter
  Double xp, yp, deltaUT, LOD, ddpsi, ddeps;
  eop(timeGPS, xp, yp, deltaUT, LOD, ddpsi, ddeps);

  // The IAU 1980 Theroy of Nutation
  Double dpsi, deps, Omega;
  nutationIAU1980(timeGPS, dpsi, deps, Omega);
  dpsi += ddpsi; // Korr. mit dpsi aus EOP
  deps += ddeps;

  const Double JC     = timeGPS2JC(timeGPS);
  const Double zetaA  = (2306.2181*JC + 0.30188*JC*JC + 0.017998*JC*JC*JC)* DEG2RAD/3600;
  const Double thetaA = (2004.3109*JC - 0.42665*JC*JC - 0.041883*JC*JC*JC)* DEG2RAD/3600;
  const Double zA     = (2306.2181*JC + 1.09468*JC*JC + 0.018203*JC*JC*JC)* DEG2RAD/3600;
  const Double epsA   = (84381.448 -46.8150*JC - 0.00059*JC*JC + 0.001813*JC*JC*JC) * DEG2RAD/3600;

  // GST
  const Double theta = GST(timeGPS2UTC(timeGPS)+seconds2time(deltaUT), dpsi, epsA, Omega);

  return PolMatrix(xp,yp) *
         RotMatrix(theta) *
         NutationMatrix(epsA,dpsi,deps) *
         PraezessionMatrix(zetaA,thetaA,zA);
}

/***********************************************/

void EarthRotationIers1996::earthOrientationParameter(const Time &timeGPS, Double &xp, Double &yp, Double &sp, Double &deltaUT, Double &LOD, Double &X, Double &Y, Double &S) const
{
  Double ddpsi, ddeps;
  eop(timeGPS, xp, yp, deltaUT, LOD, ddpsi, ddeps);
  sp = 0;
  X  = 0; // wrong: must transform from nutation1980 + dpsi,deps
  Y  = 0;
  S  = 0;
}

/***********************************************/

// Erdrotationsparameter
void EarthRotationIers1996::eop(const Time &timeGPS, Double &xp, Double &yp, Double &deltaUT, Double &LOD, Double &ddpsi, Double &ddeps) const
{
  try
  {
    // interpolate EOP file
    // --------------------
    xp = yp = deltaUT = LOD = 0;
    ddpsi = ddeps = 0;
    if(EOP.size())
    {
      const Time timeUTC = timeGPS2UTC(timeGPS);
      Matrix eop = polynomial.interpolate({timeUTC}, EOP);
      xp      = eop(0,0);
      yp      = eop(0,1);
      deltaUT = eop(0,2) + (timeGPS-timeUTC).seconds();
      LOD     = eop(0,3);
      ddpsi   = eop(0,4);
      ddeps   = eop(0,5);
    }

#ifndef GROOPS_DISABLE_IERS
    // Models
    // ------
    // diurnal and semidiurnal variations in EOP (x,y,UT1) from ocean tides
    const Double mjdUTC = timeGPS2UTC(timeGPS).mjd();
    Double corx=0, cory=0, cort=0;
    ray(mjdUTC, corx, cory, cort);
    xp      += corx * DEG2RAD/3600;
    yp      += cory * DEG2RAD/3600;
    deltaUT += cort;
#endif
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

Double EarthRotationIers1996::GMST(const Time &timeUT1) const
{
  Double Tu0 = (timeUT1.mjdInt()-51544.5)/36525.0;

  Double GMST0 = (6.0/24 + 41.0/(24*60) + 50.54841/(24*60*60))
               + (8640184.812866/(24*60*60))*Tu0
               + (0.093104/(24*60*60))*Tu0*Tu0
               + (-6.2e-6/(24*60*60))*Tu0*Tu0*Tu0;
  Double r     = 1.002737909350795 + 5.9006e-11*Tu0 - 5.9e-15*Tu0*Tu0;
  return fmod(2*PI*(GMST0 + r * timeUT1.mjdMod()), 2*PI);
}

/***********************************************/

Double EarthRotationIers1996::GST(const Time &timeUT1, Double dpsi, Double epsA, Double Omega) const
{
  return fmod(GMST(timeUT1) + dpsi*cos(epsA)+0.00264*DEG2RAD/3600*sin(Omega)+0.000063*DEG2RAD/3600*sin(2*Omega), 2*PI);
}

/***********************************************/

void EarthRotationIers1996::nutationIAU1980(const Time &timeGPS, Double &dpsi, Double &deps, Double &Omega) const
{
  const Double JC  = timeGPS2JC(timeGPS);
  const Vector F   = Planets::fundamentals(timeGPS);
  const Vector arg = argument * F;

  Omega = F(4);
  dpsi  = 0;
  deps  = 0;
  for(UInt i=0; i<argument.rows(); i++)
  {
    dpsi += (psiFactor(i,0)+psiFactor(i,1)*JC)*sin(arg(i));
    deps += (epsFactor(i,0)+epsFactor(i,1)*JC)*cos(arg(i));
  }
}

/***********************************************/
/***********************************************/

Rotary3d EarthRotationIers1996::RotMatrix(Double GST) const
{
  return rotaryZ(Angle(GST));
}

/***********************************************/

Rotary3d EarthRotationIers1996::PolMatrix(Double xp, Double yp) const
{
  return rotaryY(Angle(-xp)) * rotaryX(Angle(-yp));
}

/***********************************************/

Rotary3d EarthRotationIers1996::PraezessionMatrix(Double zetaA, Double thetaA, Double zA) const
{
  return rotaryZ(Angle(-zA))*rotaryY(Angle(thetaA))*rotaryZ(Angle(-zetaA));
}

/***********************************************/

Rotary3d EarthRotationIers1996::NutationMatrix(Double epsA, Double dpsi, Double deps) const
{
  return rotaryX(Angle(-epsA-deps))*rotaryZ(Angle(-dpsi))*rotaryX(Angle(epsA));
}

/***********************************************/
/***********************************************/
