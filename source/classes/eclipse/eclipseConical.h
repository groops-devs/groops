/***********************************************/
/**
* @file eclipseConical.h
*
* @brief Shadowing of satellites by moon and Earth.
* @see Eclipse
*
* @author Torsten Mayer-Guerr
* @date 2020-03-08
*
*/
/***********************************************/

#ifndef __GROOPS_ECLIPSECONICAL__
#define __GROOPS_ECLIPSECONICAL__

// Latex documentation
#ifdef DOCSTRING_Eclipse
static const char *docstringEclipseConical = R"(
\subsection{Conical}
\fig{!hb}{0.8}{eclipseConical}{fig:eclipseConical}{Modelling umbra and penumbra.}
)";
#endif

/***********************************************/

#include "classes/eclipse/eclipse.h"

/***** CLASS ***********************************/

/** @brief Shadowing of satellites by moon and Earth.
* @ingroup eclipseGroup
* @see Eclipse */
class EclipseConical : public Eclipse
{
public:
  EclipseConical(Config &/*config*/) {}

  virtual Double factor(const Time &timeGPS, const Vector3d &position, EphemeridesPtr ephemerides) const override;
};

/***********************************************/

inline Double EclipseConical::factor(const Time &timeGPS, const Vector3d &position, EphemeridesPtr ephemerides) const
{
  try
  {
    if(!ephemerides)
      throw(Exception("No ephemerides given"));

    const Vector3d posSun   = ephemerides->position(timeGPS, Ephemerides::SUN);
    const Vector3d posMoon  = ephemerides->position(timeGPS, Ephemerides::MOON);
    return shadowScalingFactor(position, posSun, Vector3d(), R_Earth)
         * shadowScalingFactor(position, posSun, posMoon,    R_Moon);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
