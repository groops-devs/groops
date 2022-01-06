/***********************************************/
/**
* @file earthRotationStarCamera.h
*
* @brief Earth rotation from quaternion (star camera) file.
* @see EarthRotation
*
* @author Andreas Kvas
* @date 2019-01-22
*
*/
/***********************************************/

#ifndef __GROOPS_EARTHROTATIONSTARCAMAMERA__
#define __GROOPS_EARTHROTATIONSTARCAMAMERA__

// Latex documentation
#ifdef DOCSTRING_EarthRotation
static const char *docstringEarthRotationStarCamera = R"(
\subsection{StarCamera}
This class reads quaternions from an instrument file and interpolates to the given time stamp.
)";
#endif

/***********************************************/

#include "base/polynomial.h"
#include "classes/earthRotation/earthRotation.h"

/***** CLASS ***********************************/

/** @brief Earth rotation from quaternion input.
* @ingroup earthRotationGroup
* @see EarthRotation */
class EarthRotationStarCamera : public EarthRotation
{
  Polynomial polynomial;
  Matrix     quaternions;

public:
  EarthRotationStarCamera(Config &config);

  void earthOrientationParameter(const Time &/*timeGPS*/, Double &/*xp*/, Double &/*yp*/, Double &/*sp*/, Double &/*deltaUT*/, Double &/*LOD*/, Double &/*X*/, Double &/*Y*/, Double &/*S*/) const { throw Exception("Earth orientation parameters can not be determined from quaterions alone.");}
  Rotary3d rotaryMatrix(const Time &timeGPS) const;
};

/***********************************************/

#endif /* __GROOPS_EARTHROTATION__ */
