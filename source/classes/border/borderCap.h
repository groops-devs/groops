/***********************************************/
/**
* @file borderCap.h
*
* @brief Spherical cap.
* @see Border
*
* @author Annette Eicker
* @date 2004-10-28
*
*/
/***********************************************/

#ifndef __GROOPS_BORDERCAP__
#define __GROOPS_BORDERCAP__

// Latex documentation
#ifdef DOCSTRING_Border
static const char *docstringBorderCap = R"(
\subsection{Cap}
The region is defined by a spherical cap with the center given in geographical coordinates
longitude (\config{lambdaCenter}) and latitude (\config{phiCenter}).
The radius of the cap is given as aperture angle \config{psi}.

\fig{!hb}{0.4}{borderCap}{fig:borderCap}{spherical cap}
)";
#endif

/***********************************************/

#include "config/config.h"
#include "classes/border/border.h"

/***** CLASS ***********************************/

/** @brief Spherical cap.
* @ingroup borderGroup
* @see Border */
class BorderCap : public BorderBase
{
  Double   cosPsi;
  Vector3d center;
  Bool     exclude;

public:
  BorderCap(Config &config);

  Bool isInnerPoint(Angle lambda, Angle phi) const;
  Bool isExclude() const {return exclude;}
};

/***********************************************/

inline BorderCap::BorderCap(Config &config)
{
  Angle psi;
  Angle lambdaCenter, phiCenter;

  readConfig(config, "lambdaCenter", lambdaCenter, Config::MUSTSET,  "",  "longitude of the center of the cap");
  readConfig(config, "phiCenter",    phiCenter,    Config::MUSTSET,  "",  "latitude of the center of the cap");
  readConfig(config, "psi",          psi,          Config::MUSTSET,  "",  "aperture angle (radius)");
  readConfig(config, "exclude",      exclude,      Config::DEFAULT,  "0", "dismiss points inside");
  if(isCreateSchema(config)) return;

  cosPsi = cos(psi);
  center = polar(lambdaCenter, phiCenter, 1.0);
}

/***********************************************/

inline Bool BorderCap::isInnerPoint(Angle lambda, Angle phi) const
{
  Vector3d pkt = polar(lambda, phi, 1.0);
  return cosPsi<=inner(pkt,center);
}

/***********************************************/

#endif /* __GROOPS_BORDER__ */
