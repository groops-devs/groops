/***********************************************/
/**
* @file earthRotationEra.h
*
* @brief Simple rotation about z-axis (ERA).
* @see EarthRotation
*
* @author Torsten Mayer-Guerr
* @date 2009-11-25
*
*/
/***********************************************/

#ifndef __GROOPS_EARTHROTATIONERA__
#define __GROOPS_EARTHROTATIONERA__

// Latex documentation
#ifdef DOCSTRING_EarthRotation
static const char *docstringEarthRotationEra = R"(
\subsection{Earth Rotation Angle (ERA)}
The transformation is realized as rotation about the z-axis.
The angle ist given by the Earth Rotation Angle (ERA).
)";
#endif

/***********************************************/

#include "base/planets.h"
#include "classes/earthRotation/earthRotation.h"

/***** CLASS ***********************************/

/** @brief Simple rotation about z-axis (ERA).
* @ingroup earthRotationGroup
* @see EarthRotation */
class EarthRotationEra : public EarthRotation
{
public:
  inline EarthRotationEra(Config &config);

  inline Rotary3d rotaryMatrix (const Time &timeGPS) const;
  inline Vector3d rotaryAxis   (const Time &timeGPS) const;
  inline Vector3d rotaryAxisDerivate(const Time &timeGPS) const;
};

/***********************************************/

inline EarthRotationEra::EarthRotationEra(Config &/*config*/)
{
}

/***********************************************/

inline Rotary3d EarthRotationEra::rotaryMatrix(const Time &timeGPS) const
{
  return rotaryZ(Angle(Planets::ERA(timeGPS2UTC(timeGPS))));
}

/***********************************************/

inline Vector3d EarthRotationEra::rotaryAxis(const Time &/*timeGPS*/) const
{
  return Vector3d(0.0, 0.0, 2*PI*1.00273781191135448/86400);
}

/***********************************************/

inline Vector3d EarthRotationEra::rotaryAxisDerivate(const Time &/*timeGPS*/) const
{
  return Vector3d(0.0, 0.0, 0.0);
}

/***********************************************/

#endif /* __GROOPS_EARTHROTATION__ */
