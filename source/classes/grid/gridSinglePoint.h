/***********************************************/
/**
* @file gridSinglePoint.h
*
* @brief Single point.
* @see Grid
* single point at sphere/ellipsoid.
*
* @author Torsten Mayer-Guerr
* @date 2010-11-01
*
*/
/***********************************************/

#ifndef __GROOPS_GRIDSINGLEPOINT__
#define __GROOPS_GRIDSINGLEPOINT__

// Latex documentation
#ifdef DOCSTRING_Grid
static const char *docstringGridSinglePoint = R"(
\subsection{SinglePoint}
Creates one single point.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "classes/grid/grid.h"

/***** CLASS ***********************************/

/** @brief Single point.
* @ingroup gridGroup
* @see Grid */
class GridSinglePoint : public GridBase
{
public:
  GridSinglePoint(Config &config);
};

/***********************************************/

inline GridSinglePoint::GridSinglePoint(Config &config)
{
  try
  {
    Angle     L, B;
    Double    a, f, height, area = NAN_EXPR;

    readConfig(config, "L",                 L,          Config::MUSTSET,  "",   "longitude");
    readConfig(config, "B",                 B,          Config::MUSTSET,  "",   "latitude");
    readConfig(config, "height",            height,     Config::DEFAULT,  "0.", "ellipsoidal height");
    readConfig(config, "area",              area,       Config::OPTIONAL, "",   "associated area element on unit sphere");
    readConfig(config, "R",                 a,          Config::DEFAULT,  STRING_DEFAULT_GRS80_a, "major axsis of the ellipsoid/sphere");
    readConfig(config, "inverseFlattening", f,          Config::DEFAULT,  STRING_DEFAULT_GRS80_f, "flattening of the ellipsoid, 0: sphere");
    if(isCreateSchema(config)) return;

    Ellipsoid ellipsoid(a,f);
    points.push_back(ellipsoid(L, B, height));
    if(!std::isnan(area))
      areas.push_back(area);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif /* __GROOPS__ */
