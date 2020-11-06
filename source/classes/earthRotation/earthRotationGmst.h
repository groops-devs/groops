/***********************************************/
/**
* @file earthRotationGmst.h
*
* @brief Simple rotation about z-axis (GMST).
* @see EarthRotation
*
* @author Torsten Mayer-Guerr
* @date 2005-12-01
*
*/
/***********************************************/

#ifndef __GROOPS_EARTHROTATIONGMST__
#define __GROOPS_EARTHROTATIONGMST__

// Latex documentation
#ifdef DOCSTRING_EarthRotation
static const char *docstringEarthRotationGmst = R"(
\subsection{Gmst}
The transformation is realized as rotation about the z-axis.
The angle ist given by the Greenwich Mean Siderial Time (GMST).
\begin{verbatim}
  Double Tu0 = (timeUTC.mjdInt()-51544.5)/36525.0;

  Double GMST0 = (6.0/24 + 41.0/(24*60) + 50.54841/(24*60*60))
               + (8640184.812866/(24*60*60))*Tu0
               + (0.093104/(24*60*60))*Tu0*Tu0
               + (-6.2e-6/(24*60*60))*Tu0*Tu0*Tu0;
  Double r = 1.002737909350795 + 5.9006e-11*Tu0 - 5.9e-15*Tu0*Tu0;
  GMST = fmod(2*PI*(GMST0 + r * timeUTC.mjdMod()), 2*PI);
\end{verbatim}
)";
#endif

/***********************************************/

#include "base/planets.h"
#include "classes/earthRotation/earthRotation.h"

/***** CLASS ***********************************/

/** @brief Simple rotation about z-axis (GMST).
* @ingroup earthRotationGroup
* @see EarthRotation */
class EarthRotationGmst : public EarthRotation
{
public:
  inline EarthRotationGmst(Config &config);

  inline Rotary3d rotaryMatrix (const Time &timeGPS) const;
  inline Vector3d rotaryAxis   (const Time &timeGPS) const;
  inline Vector3d rotaryAxisDerivate(const Time &timeGPS) const;
};

/***********************************************/

inline EarthRotationGmst::EarthRotationGmst(Config &/*config*/)
{
}

/***********************************************/

inline Rotary3d EarthRotationGmst::rotaryMatrix(const Time &timeGPS) const
{
  return rotaryZ(Angle(Planets::gmst(timeGPS2UTC(timeGPS))));
}

/***********************************************/

inline Vector3d EarthRotationGmst::rotaryAxis(const Time &/*timeGPS*/) const
{
  return Vector3d(0., 0., 2*PI*1.002737909350795/86400);
}

/***********************************************/

inline Vector3d EarthRotationGmst::rotaryAxisDerivate(const Time &/*timeGPS*/) const
{
  return Vector3d(0., 0., 0.);
}

/***********************************************/

#endif /* __GROOPS_EARTHROTATION__ */
