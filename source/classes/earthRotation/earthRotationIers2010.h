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
