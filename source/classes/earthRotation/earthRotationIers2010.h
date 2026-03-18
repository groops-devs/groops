/***********************************************/
/**
* @file earthRotationIers2010.h
*
* @brief According to IERS2010 conventions.
* @see EarthRotation
*
* @author Torsten Mayer-Guerr
* @date 2010-03-01
*
*/
/***********************************************/

#ifndef __GROOPS_EARTHROTATIONIERS2010__
#define __GROOPS_EARTHROTATIONIERS2010__

// Latex documentation
#ifdef DOCSTRING_EarthRotation
static const char *docstringEarthRotationIers2010 = R"(
\subsection{Iers2010}
This class realize the transformation according to the IERS2010 conventions
given by the \emph{International Earth Rotation and Reference Systems Service} (IERS).
A file with the earth orientation parameter is needed (\configFile{inputfileEOP}{earthOrientationParameter}).
)";
#endif

/***********************************************/

#include "base/polynomial.h"
#include "classes/earthRotation/earthRotation.h"

/***** CLASS ***********************************/

/** @brief According to IERS2010 conventions.
 * In addtion to those interpolated values from the input file, 
 * - diurnal and semi-diurnal variations from ocean tides are considered by the IERS routine \c ORTHO_EOP.F for
 *   x-pole, y-pole and dUT1
 * - librations are considered by the routine \c PMSDNUT2.F for x-pole and y-pole
 * - librations are considered by the routine \c UTLIBR.F for dUT1 and LOD 
 * 
 * sp is calculated with the ERFA function \c eraSp00. X,Y, and S are calculated either by the truncated 
 * precession & nutation model with the ERFA function \c eraXys00b or by the untruncated 
 * precession & nutation model with the ERFA functions \c eraXy06 and \c eraS06.
* @ingroup earthRotationGroup
* @see EarthRotation */
class EarthRotationIers2010 : public EarthRotation
{
  Bool              useTruncated;
  Polynomial        polynomial;
  Matrix            EOP;

public:
  EarthRotationIers2010(Config &config);

  void earthOrientationParameter(const Time &timeGPS, Double &xp, Double &yp, Double &sp, Double &deltaUT, Double &LOD, Double &X, Double &Y, Double &S) const;
};

/***********************************************/

#endif /* __GROOPS_EARTHROTATION__ */
