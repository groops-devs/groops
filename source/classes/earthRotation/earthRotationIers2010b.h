/***********************************************/
/**
* @file earthRotationIers2010b.h
*
* @brief According to IERS2010 conventions (with HF EOP model).
* @see EarthRotation
*
* @author Torsten Mayer-Guerr
* @date 2019-05-15
*
*/
/***********************************************/

#ifndef __GROOPS_EARTHROTATIONIERS2010b__
#define __GROOPS_EARTHROTATIONIERS2010b__

// Latex documentation
#ifdef DOCSTRING_EarthRotation
static const char *docstringEarthRotationIers2010b = R"(
\subsection{Iers2010b}\label{earthRotationType:iers2010b}
This class realize the transformation according to the IERS2010 conventions
given by the \emph{International Earth Rotation and Reference Systems Service} (IERS).
A file with the earth orientation parameter is needed (\configFile{inputfileEOP}{earthOrientationParameter}).
Includes additional high-frequency EOP models (\configFile{inputfileDoodsonEOP}{doodsonEarthOrientationParameter}).
)";
#endif

/***********************************************/

#include "base/polynomial.h"
#include "files/fileDoodsonEarthOrientationParameter.h"
#include "classes/earthRotation/earthRotation.h"

/***** CLASS ***********************************/

/** @brief According to IERS2010 conventions.
* @ingroup earthRotationGroup
* @see EarthRotation */
class EarthRotationIers2010b : public EarthRotation
{
  Polynomial        polynomial;
  Matrix            EOP;
  DoodsonEop        doodsonEop;
  Matrix            doodsonMatrix;

public:
  EarthRotationIers2010b(Config &config);

  void earthOrientationParameter(const Time &timeGPS, Double &xp, Double &yp, Double &sp, Double &deltaUT, Double &LOD, Double &X, Double &Y, Double &S) const;
};

/***********************************************/

#endif /* __GROOPS_EARTHROTATION__ */
