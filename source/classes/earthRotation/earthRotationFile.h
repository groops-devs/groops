/***********************************************/
/**
* @file earthRotationFile.h
*
* @brief Interpolated values from file.
* @see EarthRotation
*
* @author Torsten Mayer-Guerr
* @date 2010-03-01
*
*/
/***********************************************/

#ifndef __GROOPS_EARTHROTATIONFILE__
#define __GROOPS_EARTHROTATIONFILE__

// Latex documentation
#ifdef DOCSTRING_EarthRotation
static const char *docstringEarthRotationFile = R"(
\subsection{File}\label{earthRotationType:file}
This class realize the transformation by interpolation from file.
This file can be created with \program{EarthOrientationParameterTimeSeries}.
)";
#endif

/***********************************************/

#include "base/polynomial.h"
#include "classes/earthRotation/earthRotation.h"

/***** CLASS ***********************************/

/** @brief According to IERS2010 conventions.
* @ingroup earthRotationGroup
* @see EarthRotation */
class EarthRotationFile : public EarthRotation
{
  Polynomial        polynomial;
  Matrix            EOP;

public:
  EarthRotationFile(Config &config);

  void earthOrientationParameter(const Time &timeGPS, Double &xp, Double &yp, Double &sp, Double &deltaUT, Double &LOD, Double &X, Double &Y, Double &S) const;
};

/***********************************************/

#endif /* __GROOPS_EARTHROTATION__ */
