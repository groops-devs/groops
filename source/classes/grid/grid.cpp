/***********************************************/
/**
* @file grid.cpp
*
* @brief point distributions on the sphere/ellipsoid.
*
* @author Torsten Mayer-Guerr
* @date 2002-05-31
*
*/
/***********************************************/

#define DOCSTRING_Grid

#include "base/import.h"
#include "config/configRegister.h"
#include "classes/grid/gridGeograph.h"
#include "classes/grid/gridTriangleVertex.h"
#include "classes/grid/gridTriangleCenter.h"
#include "classes/grid/gridGauss.h"
#include "classes/grid/gridReuter.h"
#include "classes/grid/gridCorput.h"
#include "classes/grid/gridDriscoll.h"
#include "classes/grid/gridSinglePoint.h"
#include "classes/grid/gridSinglePointCartesian.h"
#include "classes/grid/gridFile.h"
#include "classes/grid/grid.h"

/***********************************************/

GROOPS_REGISTER_CLASS(Grid, "gridType",
                      GridGeograph,
                      GridTriangleVertex,
                      GridTriangleCenter,
                      GridGauss,
                      GridReuter,
                      GridCorput,
                      GridDriscoll,
                      GridSinglePoint,
                      GridSinglePointCartesian,
                      GridFile)

GROOPS_READCONFIG_UNBOUNDED_CLASS(Grid, "gridType")

/***********************************************/

Grid::Grid(Config &config, const std::string &name)
{
  try
  {
    std::unique_ptr<GridBase> base;

    std::string choice;
    while(readConfigChoice(config, name, choice, Config::OPTIONAL, "", "point distributions on the sphere/ellipsoid"))
    {
      if(readConfigChoiceElement(config, "geograph",       choice, "along lines of geographical coordinates"))
        base = std::unique_ptr<GridBase>(new GridGeograph(config));
      if(readConfigChoiceElement(config, "triangleVertex", choice, "triangle grid (vertcies)"))
        base = std::unique_ptr<GridBase>(new GridTriangleVertex(config));
      if(readConfigChoiceElement(config, "triangleCenter", choice, "triangle grid (center points)"))
        base = std::unique_ptr<GridBase>(new GridTriangleCenter(config));
      if(readConfigChoiceElement(config, "gauss",          choice, "gauss-legendre grid"))
        base = std::unique_ptr<GridBase>(new GridGauss(config));
      if(readConfigChoiceElement(config, "reuter",         choice, "reuter grid"))
        base = std::unique_ptr<GridBase>(new GridReuter(config));
      if(readConfigChoiceElement(config, "corput",         choice, "pseudo random distribution"))
        base = std::unique_ptr<GridBase>(new GridCorput(config));
      if(readConfigChoiceElement(config, "driscoll",       choice, "driscoll-healy grid"))
        base = std::unique_ptr<GridBase>(new GridDriscoll(config));
      if(readConfigChoiceElement(config, "singlePoint",    choice, "one single point"))
        base = std::unique_ptr<GridBase>(new GridSinglePoint(config));
      if(readConfigChoiceElement(config, "singlePointCartesian", choice, "one single point"))
        base = std::unique_ptr<GridBase>(new GridSinglePointCartesian(config));
      if(readConfigChoiceElement(config, "file",           choice, "read points (with values) from file"))
        base = std::unique_ptr<GridBase>(new GridFile(config));
      endChoice(config);
      if(isCreateSchema(config))
        return;

      points_.insert(points_.end(), base->points.begin(), base->points.end());
      areas_.insert (areas_.end(),  base->areas.begin(),  base->areas.end());
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
