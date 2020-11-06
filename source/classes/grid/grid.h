/***********************************************/
/**
* @file grid.h
*
* @brief Point distributions on the sphere/ellipsoid.
*
* @author Torsten Mayer-Guerr
* @date 2002-05-31
*
*/
/***********************************************/

#ifndef __GROOPS_GRID__
#define __GROOPS_GRID__

// Latex documentation
#ifdef DOCSTRING_Grid
static const char *docstringGrid = R"(
\section{Grid}\label{gridType}
This class generates a set of grid points. In a first step, the grid
is always generated globally, with \configClass{border}{borderType} a regional
subset of points can be extracted from the global grid. The parameters
\config{R} and \config{inverseFlattening} define the shape of the ellipsoid
on which the grid is generated. In case \config{inverseFlattening} is
chosen as zero, a sphere is used. With \config{height} the distance of
the points above the ellipsoid can be defined. In addition to the location
of the points, weights are assigned to each of the points. These weights
can be regarded as the surface element associated with each grid point.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"

/**
* @defgroup gridGroup Grid
* @brief Point distributions on the sphere/ellipsoid.
* @ingroup classesGroup
* The interface is given by @ref Grid.
* An Instance can be created by @ref readConfig. */
/// @{

/***** TYPES ***********************************/

class Grid;
typedef std::shared_ptr<Grid> GridPtr;

/***** CLASS ***********************************/

/** @brief Point distributions on the sphere/ellipsoid.
* An instance of this class can be created with @ref readConfig. */
class Grid
{
  std::vector<Vector3d> points_;
  std::vector<Double>   areas_;

public:
  /// Constructor.
  Grid(Config &config, const std::string &name);

  /** @brief Point distribution. */
  const std::vector<Vector3d> &points() const {return points_;}

  /** @brief Area element relating to each point.
  * Computed on a unit sphere (only approximative). */
  const std::vector<Double> &areas() const {return areas_;}

  /** @brief creates an derived instance of this class. */
  static GridPtr create(Config &config, const std::string &name) {return GridPtr(new Grid(config, name));}
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class Grid.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and an class without points is created.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] grid Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates Grid */
template<> Bool readConfig(Config &config, const std::string &name, GridPtr &grid, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***** CLASS ***********************************/

// Internal class
class GridBase
{
public:
  std::vector<Vector3d> points;
  std::vector<Double>   areas;

  virtual ~GridBase() {}
};

/***********************************************/

#endif /* __GROOPS_GRID__ */

