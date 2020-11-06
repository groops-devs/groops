/***********************************************/
/**
* @file gridSinglePointCartesian.h
*
* @brief Single point.
* @see Grid
* single point.
*
* @author Torsten Mayer-Guerr
* @date 2014-09-25
*
*/
/***********************************************/

#ifndef __GROOPS_GRIDSINGLEPOINTCARTESIAN__
#define __GROOPS_GRIDSINGLEPOINTCARTESIAN__

// Latex documentation
#ifdef DOCSTRING_Grid
static const char *docstringGridSinglePointCartesian = R"(
\subsection{SinglePointCartesian}
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
class GridSinglePointCartesian : public GridBase
{
public:
  GridSinglePointCartesian(Config &config);
};

/***********************************************/

inline GridSinglePointCartesian::GridSinglePointCartesian(Config &config)
{
  try
  {
    Double x, y, z, area = NAN_EXPR;

    readConfig(config, "x",    x,    Config::MUSTSET,  "", "[m]");
    readConfig(config, "y",    y,    Config::MUSTSET,  "", "[m]");
    readConfig(config, "z",    z,    Config::MUSTSET,  "", "[m]");
    readConfig(config, "area", area, Config::OPTIONAL, "", "associated area element on unit sphere");
    if(isCreateSchema(config)) return;

    points.push_back(Vector3d(x,y,z));
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
