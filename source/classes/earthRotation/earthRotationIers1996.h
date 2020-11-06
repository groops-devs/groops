/***********************************************/
/**
* @file earthRotationIers1996.h
*
* @brief According to IERS1996 conventions.
* @see EarthRotation
*
* @author Torsten Mayer-Guerr
* @date 2003-09-20
*
*/
/***********************************************/

#ifndef __GROOPS_EARTHROTATIONIERS1996__
#define __GROOPS_EARTHROTATIONIERS1996__

// Latex documentation
#ifdef DOCSTRING_EarthRotation
static const char *docstringEarthRotationIers1996 = R"(
\subsection{Iers1996}
Very old.
)";
#endif

/***********************************************/

#include "base/polynomial.h"
#include "classes/earthRotation/earthRotation.h"

/***** CLASS ***********************************/

/** @brief According to IERS1996 conventions.
* @ingroup earthRotationGroup
* @see EarthRotation */
class EarthRotationIers1996 : public EarthRotation
{
  Polynomial        polynomial;
  std::vector<Time> times; // UTC
  Matrix            EOP;
  Matrix            argument, psiFactor, epsFactor; // nutation series

  void     eop(const Time &timeGPS, Double &xp, Double &yp, Double &deltaUT, Double &LOD, Double &ddpsi, Double &ddeps) const;
  Double   GST (const Time &timeUT1, Double dpsi, Double epsA, Double Omega) const;
  Double   GMST(const Time &timeUT1)                                         const;
  void     nutationIAU1980(const Time &timeGPS, Double &dpsi, Double &deps, Double &Omega) const;
  Rotary3d RotMatrix(Double GST) const;
  Rotary3d PolMatrix(Double xp, Double yp) const;
  Rotary3d PraezessionMatrix(Double zetaA, Double thetaA, Double zA) const;
  Rotary3d NutationMatrix(Double epsA, Double dpsi, Double deps) const;

public:
  EarthRotationIers1996(Config &config);

  Rotary3d rotaryMatrix(const Time &timeGPS) const;
  void     earthOrientationParameter(const Time &timeGPS, Double &xp, Double &yp, Double &sp, Double &deltaUT, Double &LOD, Double &X, Double &Y, Double &S) const;
};

/***********************************************/

#endif /* __GROOPS_EARTHROTATION__ */
