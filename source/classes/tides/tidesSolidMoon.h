/***********************************************/
/**
* @file tidesSolidMoon.h
*
* @brief Solid moon tides.
* @see Tides
*
* @author Beate Klinger
* @date 2013-xx-xx
*
*/
/***********************************************/


#ifndef __GROOPS_TIDESSOLIDMOON__
#define __GROOPS_TIDESSOLIDMOON__

// Latex documentation
#ifdef DOCSTRING_Tides
static const char *docstringTidesSolidMoon = R"(
\subsection{SolidMoonTide}
This class computes the solid moon tide according to the IERS2010 conventions.
The values of solid Moon tide external potential Love numbers are given and
there are no frequency dependent corrections of these values.
The computed result is multiplied with \config{factor}.
)";
#endif

/***********************************************/

#include "base/matrix.h"
#include "classes/tides/tides.h"

/***** CLASS ***********************************/

/** @brief Solid moon tides.
* @ingroup tidesGroup
* @see Tides */
class TidesSolidMoon : public TidesBase
{
  Double kReal20;
  Double kReal30;
  Double factor;

  void moonCoefficients1(Double GM_third, const Vector3d &third, Matrix &cnm, Matrix &snm) const;

public:
  TidesSolidMoon(Config &config);

  SphericalHarmonics sphericalHarmonics(const Time &time, const Rotary3d &rotEarth, EarthRotationPtr rotation, EphemeridesPtr ephemerides,
                                        UInt maxDegree, UInt minDegree, Double GM, Double R) const override;
};

/***********************************************/

#endif
