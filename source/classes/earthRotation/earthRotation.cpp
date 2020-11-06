/***********************************************/
/**
* @file earthRotation.cpp
*
* @brief Transformation between CRF and TRF.
*
* @author Torsten Mayer-Guerr
* @date 2001-11-12
*
*/
/***********************************************/

#define DOCSTRING_EarthRotation

#include "base/import.h"
#include "base/planets.h"
#include "config/configRegister.h"
#include "classes/earthRotation/earthRotationFile.h"
#include "classes/earthRotation/earthRotationIers2010.h"
#include "classes/earthRotation/earthRotationIers2010b.h"
#include "classes/earthRotation/earthRotationIers2003.h"
#include "classes/earthRotation/earthRotationIers1996.h"
#include "classes/earthRotation/earthRotationGmst.h"
#include "classes/earthRotation/earthRotationEra.h"
#include "classes/earthRotation/earthRotationZAxis.h"
#include "classes/earthRotation/moonRotation.h"
#include "classes/earthRotation/earthRotationStarCamera.h"
#include "classes/earthRotation/earthRotation.h"


/***********************************************/

GROOPS_REGISTER_CLASS(EarthRotation, "earthRotationType",
                      EarthRotationFile,
                      EarthRotationIers2010,
                      EarthRotationIers2010b,
                      EarthRotationIers2003,
                      EarthRotationIers1996,
                      EarthRotationGmst,
                      EarthRotationEra,
                      EarthRotationZAxis,
                      EarthRotationStarCamera,
                      MoonRotation)

GROOPS_READCONFIG_CLASS(EarthRotation, "earthRotationType")

/***********************************************/

EarthRotationPtr EarthRotation::create(Config &config, const std::string &name)
{
  try
  {
    EarthRotationPtr earthRotation;
    std::string choice;
    readConfigChoice(config, name, choice, Config::MUSTSET, "", "transformation from CRF to TRF");

    renameDeprecatedChoice(config, choice, "itrf2010",  "iers2010",  date2time(2020, 8, 24));
    renameDeprecatedChoice(config, choice, "itrf2010b", "iers2010b", date2time(2020, 8, 24));
    renameDeprecatedChoice(config, choice, "itrf2003",  "iers2003",  date2time(2020, 8, 24));

    if(readConfigChoiceElement(config, "file",     choice, "interpolated values from file"))
      earthRotation = EarthRotationPtr(new EarthRotationFile(config));
    if(readConfigChoiceElement(config, "iers2010",  choice, "IERS conventions 2010"))
      earthRotation = EarthRotationPtr(new EarthRotationIers2010(config));
    if(readConfigChoiceElement(config, "iers2010b", choice, "IERS conventions 2010 with HF EOP model"))
      earthRotation = EarthRotationPtr(new EarthRotationIers2010b(config));
    if(readConfigChoiceElement(config, "iers2003",  choice, "IERS conventions 2003"))
      earthRotation = EarthRotationPtr(new EarthRotationIers2003(config));
    if(readConfigChoiceElement(config, "itrf1996", choice, "IERS conventions 1996"))
      earthRotation = EarthRotationPtr(new EarthRotationIers1996(config));
    if(readConfigChoiceElement(config, "gmst",     choice, "z-Axis rotation with gmst"))
      earthRotation = EarthRotationPtr(new EarthRotationGmst(config));
    if(readConfigChoiceElement(config, "era",     choice, "z-Axis rotation with earth rotation angle (ERA)"))
      earthRotation = EarthRotationPtr(new EarthRotationEra(config));
    if(readConfigChoiceElement(config, "zAxis",    choice, "z-Axis rotation"))
      earthRotation = EarthRotationPtr(new EarthRotationZAxis(config));
    if(readConfigChoiceElement(config, "starCamera", choice, "earth rotation from instrument file"))
      earthRotation = EarthRotationPtr(new EarthRotationStarCamera(config));
    if(readConfigChoiceElement(config, "moonRotation", choice, "Libration angles (Moon)"))
      earthRotation = EarthRotationPtr(new MoonRotation(config));
    endChoice(config);

    return earthRotation;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Rotary3d EarthRotation::rotaryMatrix(const Time &timeGPS) const
{
  try
  {
    Double xp, yp, sp, deltaUT, LOD, X, Y, S;
    earthOrientationParameter(timeGPS, xp, yp, sp, deltaUT, LOD, X, Y, S);

    const Double ERA = Planets::ERA(timeGPS2UTC(timeGPS) + seconds2time(deltaUT));
    const Double r2  = X*X + Y*Y;
    const Double E   = (r2 != 0.) ? std::atan2(Y, X) : 0.;
    const Double D   = std::atan(std::sqrt(r2/(1-r2)));

    return  rotaryX(Angle(-yp)) * rotaryY(Angle(-xp)) *
            rotaryZ(Angle(sp+ERA-S-E)) *
            rotaryY(Angle(D)) * rotaryZ(Angle(E));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector3d EarthRotation::rotaryAxis(const Time &timeGPS) const
{
  try
  {
    const UInt   degree = 8;
    const Double dt     = 30*60; // 30 min
    std::vector<Rotary3d> rot(degree+1);
    for(UInt i=0; i<=degree; i++)
      rot.at(i) = rotaryMatrix(timeGPS + seconds2time((i-degree/2.)*dt));

    Rotary3d rot0 = inverse(rot.at(degree/2));
    for(UInt i=0; i<=degree; i++)
      rot.at(i) = rot0*rot.at(i);

    std::vector<Double> coeff = {1./280, -4./105, 1./5., -4./5., 0, 4./5., -1./5., 4./105., -1./280};
    Vector3d omega;
    for(UInt i=0; i<coeff.size(); i++)
    {
      const Matrix R = rot.at(i).matrix();
      omega.x() += coeff[i] * (R(1,2)-R(2,1));
      omega.y() += coeff[i] * (R(2,0)-R(0,2));
      omega.z() += coeff[i] * (R(0,1)-R(1,0));
    }

    return 0.5/dt * omega;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector3d EarthRotation::rotaryAxisDerivate(const Time &timeGPS) const
{
  const Double dt = 30*60; // 30 min
  return 0.5/dt * (rotaryAxis(timeGPS+seconds2time(dt))-rotaryAxis(timeGPS-seconds2time(dt)));
}

/***********************************************/

void EarthRotation::earthOrientationParameter(const Time &/*timeGPS*/, Double &xp, Double &yp, Double &sp, Double &deltaUT, Double &LOD, Double &X, Double &Y, Double &S) const
{
  xp = yp = sp = deltaUT = LOD = X = Y = S = 0.0;
}

/***********************************************/
