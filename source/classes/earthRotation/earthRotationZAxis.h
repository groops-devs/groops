/***********************************************/
/**
* @file earthRotationZAxis.h
*
* @brief Simple rotation about z-axis (GMST).
* @see EarthRotation
*
* @author Torsten Mayer-Guerr
* @date 2001-11-12
*
*/
/***********************************************/

#ifndef __GROOPS_EARTHROTATIONZAXIS__
#define __GROOPS_EARTHROTATIONZAXIS__

// Latex documentation
#ifdef DOCSTRING_EarthRotation
static const char *docstringEarthRotationZAxis = R"(
\subsection{Z-Axis}
The transformation is realized as rotation about the z-axis.
You must specify the angle (\config{initialAngle}) at \config{time0} and
the angular velocity (\config{angularVelocity}).
)";
#endif

/***********************************************/

#include "classes/earthRotation/earthRotation.h"

/***** CLASS ***********************************/

/** @brief Simple rotation about z-axis.
* @ingroup earthRotationGroup
* @see EarthRotation */
class EarthRotationZAxis : public EarthRotation
{
  Double    angle0;           // Drehwinkel zum Zeitpunkt 0 [rad]
  Double    angleVelocity;    // Drehung pro Sekunde [rad/sek]
  Time      time0;

public:
  inline EarthRotationZAxis(Config &config);

  inline Rotary3d rotaryMatrix (const Time &timeGPS) const;
  inline Vector3d rotaryAxis   (const Time &timeGPS) const;
  inline Vector3d rotaryAxisDerivate(const Time &timeGPS) const;
};

/***********************************************/

inline EarthRotationZAxis::EarthRotationZAxis(Config &config)
{
  readConfig(config, "initialAngle",    angle0,        Config::MUSTSET, "5.133658456", "Angle at time0 [rad]");
  readConfig(config, "angularVelocity", angleVelocity, Config::MUSTSET, "7.29211585531e-5", "[rad/s]");
  readConfig(config, "time0",           time0,         Config::MUSTSET, "51740.5", "");
  if(isCreateSchema(config)) return;
}

/***********************************************/

inline Rotary3d EarthRotationZAxis::rotaryMatrix(const Time &timeGPS) const
{
  return rotaryZ(Angle(angle0 + angleVelocity*(timeGPS-time0).seconds()));
}

/***********************************************/

inline Vector3d EarthRotationZAxis::rotaryAxis(const Time &/*timeGPS*/) const
{
  return Vector3d(0.0, 0.0, angleVelocity);
}

/***********************************************/

inline Vector3d EarthRotationZAxis::rotaryAxisDerivate(const Time &/*timeGPS*/) const
{
  return Vector3d(0.0, 0.0, 0.0);
}

/***********************************************/

#endif /* __GROOPS_EARTHROTATION__ */
