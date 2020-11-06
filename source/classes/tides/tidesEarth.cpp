/***********************************************/
/**
* @file tidesEarth.cpp
*
* @brief Earth tides.
* Following the IERS conventions.
* @see Tides
*
* @author Torsten Mayer-Guerr
* @date 2002-12-13
*
*/
/***********************************************/

#include "base/import.h"
#include "base/sphericalHarmonics.h"
#include "base/doodson.h"
#include "config/config.h"
#include "files/fileEarthTide.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/tides/tidesEarth.h"

/***********************************************/

TidesEarth::TidesEarth(Config &config)
{
  try
  {
    FileName earthTidesName;
    readConfig(config, "inputfileEarthtide",   earthTidesName,       Config::MUSTSET, "{groopsDataDir}/tides/earthAnelastic2003.xml", "");
    readConfig(config, "includePermanentTide", includePermanentTide, Config::DEFAULT, "0",   "results in FALSE: zero tide, TRUE: tide free gravity field");
    readConfig(config, "factor",               factor,               Config::DEFAULT, "1.0", "the result is multplied by this factor, set -1 to substract the field");
    if(isCreateSchema(config)) return;

    readFileEarthTide(earthTidesName, kReal, kImag, kPlus, doodson20, doodson21, doodson22,
                      ampIp20, ampOp20, ampIp21, ampOp21, amp22,
                      h2_0, h2_2, l2_0, l2_2, l21_1, l22_1, h21_imag, l21_imag, h22_imag, l22_imag, h3, l3,
                      deformationArg21, deformationArg20,
                      dR21_ip, dR21_op, dR20_ip, dR20_op, dT21_ip, dT21_op, dT20_ip, dT20_op);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/************************************************/

// Geopotential (IERS-Convetions 1996, S.40) Step 1
// ------------------------------------------------
void TidesEarth::earthCoefficients1(Double GM_third, const Vector3d &third, Matrix &cnm, Matrix &snm) const
{
  // Kugelflaechenfunktionen berechnen
  Matrix Cnm, Snm;
  SphericalHarmonics::CnmSnm(1./R_Earth * third, 3, Cnm, Snm);

  Double factor = GM_third/GM_Earth;

  // Formel (1)
  for(UInt n=2; n<=3; n++)
    for(UInt m=0; m<=n; m++)
    {
      cnm(n,m) += factor/(2.*n+1.) * (kReal(n,m) * Cnm(n,m) + kImag(n,m) * Snm(n,m));
      snm(n,m) += factor/(2.*n+1.) * (kReal(n,m) * Snm(n,m) - kImag(n,m) * Cnm(n,m));
    }

  // Formel (4)
  for(UInt m=0; m<=2; m++)
  {
    cnm(4,m) += kPlus(2,m)/5. * factor * Cnm(2,m);
    snm(4,m) += kPlus(2,m)/5. * factor * Snm(2,m);
  }
}

/***********************************************/

// Geopotential (IERS-Convetions 1996, S.40) Step 2
// ------------------------------------------------
void TidesEarth::earthCoefficients2(const Time &time, Matrix &cnm, Matrix &snm) const
{
  Vector d = Doodson::arguments(time);

  // Korrektion fuer c20
  Vector thetaf = doodson20 * d;
  for(UInt i=0; i<thetaf.rows(); i++)
    cnm(2,0) += 1e-12 * (ampIp20(i) * cos(thetaf(i)) - ampOp20(i) * sin(thetaf(i)));

  // Korrektion fur c21 und s21
  thetaf = doodson21 * d;
  for(UInt i=0; i<thetaf.rows(); i++)
  {
    cnm(2,1) += 1e-12 * (ampIp21(i) * sin(thetaf(i)) + ampOp21(i) * cos(thetaf(i)));
    snm(2,1) += 1e-12 * (ampIp21(i) * cos(thetaf(i)) - ampOp21(i) * sin(thetaf(i)));
  }

  // Korrektion fur c22 und s22
  thetaf = doodson22 * d;
  for(UInt i=0; i<thetaf.rows(); i++)
  {
    cnm(2,2) +=  1e-12 * amp22(i) * cos(thetaf(i));
    snm(2,2) += -1e-12 * amp22(i) * sin(thetaf(i));
  }
}

/***********************************************/

SphericalHarmonics TidesEarth::sphericalHarmonics(const Time &time, const Rotary3d &rotEarth, EarthRotationPtr /*rotation*/, EphemeridesPtr ephemerides, UInt maxDegree, UInt minDegree, Double GM, Double R) const
{
  try
  {
    if(!ephemerides)
      throw(Exception("No ephemerides given"));

    Matrix cnm(5, Matrix::TRIANGULAR, Matrix::LOWER);
    Matrix snm(5, Matrix::TRIANGULAR, Matrix::LOWER);

    const Vector3d moon = rotEarth.rotate(ephemerides->position(time, Ephemerides::MOON));
    const Vector3d sun  = rotEarth.rotate(ephemerides->position(time, Ephemerides::SUN));

    earthCoefficients1(GM_Sun,  sun,  cnm, snm);
    earthCoefficients1(GM_Moon, moon, cnm, snm);
    earthCoefficients2(time, cnm, snm);

    if(!includePermanentTide)
      cnm(2,0) -= 4.4228e-8*(-0.31460)*kReal(2,0);  // Abzug Permanentgezeiten

    return SphericalHarmonics(GM_Earth, R_Earth, factor*cnm, factor*snm).get(maxDegree, minDegree, GM, R);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

Vector3d TidesEarth::deformation(const Time &time, const Vector3d &point, const Rotary3d &rotEarth, EarthRotationPtr /*rotation*/, EphemeridesPtr ephemerides,
                                 Double /*gravity*/, const Vector &/*hn*/, const Vector &/*ln*/) const
{
  try
  {
    if(!ephemerides)
      throw(Exception("No ephemerides given"));

    Vector3d displacement; // the result

    // local coordinate system
    Double   lambda = point.lambda();
    Double   phi    = point.phi();
    Vector3d up     = normalize(point);
    Vector3d east   = normalize(Vector3d(-up.y(), up.x(), 0.0));
    Vector3d north  = crossProduct(up, east);

    const Vector3d moon = rotEarth.rotate(ephemerides->position(time, Ephemerides::MOON));
    const Vector3d sun  = rotEarth.rotate(ephemerides->position(time, Ephemerides::SUN));

    // in-phase
    // --------
    Double h2 = h2_0 + h2_2*0.5*(3.*pow(sin(phi),2)-1.);
    Double l2 = l2_0 + l2_2*0.5*(3.*pow(sin(phi),2)-1.);
    displacement += deformationInPhase(GM_Sun,  sun,  point, h2, l2, h3, l3);
    displacement += deformationInPhase(GM_Moon, moon, point, h2, l2, h3, l3);

    // out-phase
    // ---------
    Double dUp    = 0;
    Double dEast  = 0;
    Double dNorth = 0;
    deformationOutPhase(GM_Sun,  sun,  lambda, phi, dUp, dEast, dNorth);
    deformationOutPhase(GM_Moon, moon, lambda, phi, dUp, dEast, dNorth);

    // frequency dependent correction
    // ------------------------------
    Vector d = Doodson::arguments(time);

    // diurnal band, equation (16)
    Vector thetaf = deformationArg21 * d;
    Double dUp21    = 0;
    Double dEast21  = 0;
    Double dNorth21 = 0;
    for(UInt i=0; i<thetaf.rows(); i++)
    {
      Double cosf = cos(thetaf(i)+lambda);
      Double sinf = sin(thetaf(i)+lambda);
      dUp21    += dR21_ip(i)*sinf + dR21_op(i)*cosf;
      dEast21  += dT21_ip(i)*cosf - dT21_op(i)*sinf;
      dNorth21 += dT21_ip(i)*sinf + dT21_op(i)*cosf;
    }
    dUp21    *= sin(2*phi);
    dEast21  *= sin(phi);
    dNorth21 *= cos(2*phi);

    // long periodic band, equation (17)
    thetaf = deformationArg20 * d;
    Double dUp20    = 0;
    Double dNorth20 = 0;
    for(UInt i=0; i<thetaf.rows(); i++)
    {
      Double cosf = cos(thetaf(i));
      Double sinf = sin(thetaf(i));
      dUp20    += dR20_ip(i)*cosf + dR20_op(i)*sinf;
      dNorth20 += dT20_ip(i)*cosf + dT20_op(i)*sinf;
    }
    dUp20    *= 1.5*pow(sin(phi),2)-0.5;
    dNorth20 *= sin(2*phi);

    displacement += (dUp+dUp21+dUp20)*up + (dEast+dEast21)*east + (dNorth+dNorth21+dNorth20)*north;

    return this->factor*displacement;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void TidesEarth::deformation(const std::vector<Time> &time, const std::vector<Vector3d> &point, const std::vector<Rotary3d> &rotEarth,
                             EarthRotationPtr rotation, EphemeridesPtr ephemerides, const std::vector<Double> &gravity, const Vector &hn, const Vector &ln,
                             std::vector<std::vector<Vector3d>> &disp) const
{
  try
  {
    for(UInt i=0; i<time.size(); i++)
      for(UInt k=0; k<point.size(); k++)
        disp.at(k).at(i) += deformation(time.at(i), point.at(k), rotEarth.at(i), rotation, ephemerides, gravity.at(k), hn, ln);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/************************************************/

Vector3d TidesEarth::deformationInPhase(Double GM_third, const Vector3d &third, const Vector3d &point, Double h2, Double l2, Double h3, Double l3) const
{
  Vector3d er = normalize(point);
  Vector3d eR = third;
  Double   R  = eR.normalize();
  Double   rR = inner(eR,er);
  Double   factor = GM_third/GM_Earth * R_Earth * pow(R_Earth/R,3);

  // displacement due to degree 2 tides, (eq.9)
  Vector3d displacement = factor * (h2*(1.5*rR*rR-0.5) * er + 3*l2*rR*(eR - rR*er));

  // displacement due to degree 3 tides, (eq.10)
  factor *= R_Earth/R;
  displacement += factor * (h3*(2.5*pow(rR,3)-1.5*rR) * er + l3*(7.5*rR*rR-1.5)*(eR - rR*er));

  return displacement;
}

/************************************************/

void TidesEarth::deformationOutPhase(Double GM_third, const Vector3d &third, Double lambda, Double phi, Double &dUp, Double &dEast, Double &dNorth) const
{
  Matrix Pnm     = SphericalHarmonics::Pnm(third.theta(), 1.0, 2);
  Double phij    = third.phi();
  Double dlambda = lambda-third.lambda();
  Double factor  = GM_third/GM_Earth * R_Earth * pow(R_Earth/third.r(),3);

  // equation (12) and (13)
  dNorth +=   -l21_1 * sin(phi) * sin(phi)  * factor * Pnm(2,1)*sqrt(3./5.)  * cos(dlambda)
         - 0.5*l22_1 * sin(phi) * cos(phi)  * factor * Pnm(2,2)*sqrt(12./5.) * cos(2*dlambda);
  dEast  +=    l21_1 * sin(phi) * cos(2*phi)* factor * Pnm(2,1)*sqrt(3./5.)  * sin(dlambda)
         - 0.5*l22_1 * sin(phi) * cos(phi)  * factor * Pnm(2,2)*sqrt(12./5.) * sin(phi) * sin(2*dlambda);
  // equation (14a) and (15a)
  dUp    += -0.75 * h21_imag * factor * sin(2*phij)*sin(2*phi)    * sin(dlambda)
         +  -0.75 * h22_imag * factor * pow(cos(phij)*cos(phi),2) * sin(2*dlambda);
  // equation (14b) and (15b)
  dNorth +=  -1.5 * l21_imag * factor * sin(2*phij)*cos(2*phi)      * sin(dlambda)
         +   0.75 * l22_imag * factor * pow(cos(phij),2)*sin(2*phi) * sin(2*dlambda);
  dEast  +=  -1.5 * l21_imag * factor * sin(2*phij)*sin(phi)        * cos(dlambda)
         +  -0.75 * l22_imag * factor * pow(cos(phij),2)*2*cos(phi) * cos(2*dlambda);
}

/***********************************************/
