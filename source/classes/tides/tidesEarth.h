/***********************************************/
/**
* @file tidesEarth.h
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

#ifndef __GROOPS_TIDESEARTH__
#define __GROOPS_TIDESEARTH__

// Latex documentation
#ifdef DOCSTRING_Tides
static const char *docstringTidesEarth = R"(
\subsection{EarthTide}\label{tidesType:earthTide}
This class computes the earth tide according to the IERS2003 conventions.
The values of solid Earth tide external potential Love numbers and
the frequency dependent corrections of these values are given in the file
\configFile{inputfileEarthtide}{earthTide}. The effect of the permanent tide is removed if
\config{includePermanentTide} is set to false.

The computed result is multiplied with \config{factor}.
)";
#endif

/***********************************************/

#include "base/matrix.h"
#include "classes/tides/tides.h"

/***** CLASS ***********************************/

/** @brief Earth tides.
* @ingroup tidesGroup
* Following the IERS conventions.
* @see Tides */
class TidesEarth : public TidesBase
{
  Bool             includePermanentTide;
  Matrix           kReal, kImag, kPlus;
  Matrix           doodson20, doodson21, doodson22;
  Vector           ampIp20, ampOp20;
  Vector           ampIp21, ampOp21;
  Vector           amp22;
  Double           factor;

  void  earthCoefficients1(Double GM_third, const Vector3d &third, Matrix &cnm, Matrix &snm) const;
  void  earthCoefficients2(const Time &time, Matrix &cnm, Matrix &snm) const;

  // deformation
  Double h2_0, h2_2;
  Double l2_0, l2_2;
  Double l21_1;
  Double l22_1;
  Double h21_imag, l21_imag;
  Double h22_imag, l22_imag;
  Double h3, l3;
  Matrix deformationArg21, deformationArg20;
  Vector dR21_ip, dR21_op;
  Vector dR20_ip, dR20_op;
  Vector dT21_ip, dT21_op;
  Vector dT20_ip, dT20_op;

  Vector3d deformationInPhase (Double GM_third, const Vector3d &third, const Vector3d &point, Double h2, Double l2, Double h3, Double l3) const;
  void     deformationOutPhase(Double GM_third, const Vector3d &third, Double lambda, Double phi, Double &dUp, Double &dEast, Double &dNorth) const;


public:
  TidesEarth(Config &config);

  SphericalHarmonics sphericalHarmonics(const Time &time, const Rotary3d &rotEarth, EarthRotationPtr rotation, EphemeridesPtr ephemerides,
                                        UInt maxDegree, UInt minDegree, Double GM, Double R) const override;
  Vector3d deformation(const Time &time, const Vector3d &point, const Rotary3d &rotEarth, EarthRotationPtr rotation, EphemeridesPtr ephemerides,
                       Double gravity, const Vector &hn, const Vector &ln) const override;
  void     deformation(const std::vector<Time> &time, const std::vector<Vector3d> &point, const std::vector<Rotary3d> &rotEarth,
                       EarthRotationPtr rotation, EphemeridesPtr ephemerides, const std::vector<Double> &gravity, const Vector &hn, const Vector &ln,
                       std::vector<std::vector<Vector3d>> &disp) const override;
};

/***********************************************/

#endif
