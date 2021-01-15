/***********************************************/
/**
* @file gravityfield2Deflections.cpp
*
* @brief Deflections of the vertical on a grid.
*
* @author Christian Pock
* @date 2012-06-05
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program computes the deflections of the vertical $\xi$ in north direction
and $\eta$ in east direction in radian of the \configClass{gravityfield}{gravityfieldType}
vector relative to the ellipsoidal normal.
The \configClass{gravityfield}{gravityfieldType} must provide the full gravity vector
inclusive the centrifugal part, see \configClass{gravityfield:tides:centrifugal}{tidesType:centrifugal}.

The values will be saved together with points expressed as ellipsoidal coordinates (longitude, latitude, height)
based on a reference ellipsoid with parameters \config{R} and \config{inverseFlattening}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/planets.h"
#include "files/fileGriddedData.h"
#include "classes/grid/grid.h"
#include "classes/gravityfield/gravityfield.h"
#include "misc/miscGriddedData.h"

/***** CLASS ***********************************/

/** @brief Deflections of the vertical on a grid.
* @ingroup programsGroup */
class Gravityfield2Deflections
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Gravityfield2Deflections, PARALLEL, "Deflections of the vertical - based on gravity field", Gravityfield)

/***********************************************/

void Gravityfield2Deflections::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName          fileNameGrid;
    GridPtr           grid;
    GravityfieldPtr   gravityfield;
    Time              time;
    Double            a, f;

    readConfig(config, "outputfileGriddedData", fileNameGrid,    Config::MUSTSET,  "", "xi (north), eta (east) [rad]");
    readConfig(config, "grid",                  grid,            Config::MUSTSET,  "", "");
    readConfig(config, "gravityfield",          gravityfield,    Config::MUSTSET,  "", "");
    readConfig(config, "time",                  time,            Config::OPTIONAL, "", "at this time the gravity field will be evaluated");
    readConfig(config, "R",                     a,               Config::DEFAULT,  STRING_DEFAULT_GRS80_a, "reference radius for ellipsoidal coordinates on output");
    readConfig(config, "inverseFlattening",     f,               Config::DEFAULT,  STRING_DEFAULT_GRS80_f, "reference flattening for ellipsoidal coordinates on output, 0: spherical coordinates");
    if(isCreateSchema(config)) return;

    // create grid
    // -----------
    Ellipsoid ellipsoid(a,f);
    std::vector<Vector3d> points = grid->points();
    std::vector<Double>   areas  = grid->areas();

    // Compute deflections
    // ------------------
    logStatus<<"compute deflections of the vertical"<<Log::endl;
    std::vector<Vector3d> g(points.size());
    Parallel::forEach(g, [&](UInt i)
    {
      return normalize(localNorthEastUp(points.at(i), ellipsoid).inverseTransform(gravityfield->gravity(time, points.at(i))));
    }, comm);

    if(Parallel::isMaster(comm))
    {
      std::vector<std::vector<Double>> field(2);
      for(UInt i=0; i<points.size(); i++)
      {
        field.at(0).at(i) = g.at(i).x();
        field.at(1).at(i) = g.at(i).y();
      }

      logStatus<<"save values <"<<fileNameGrid<<">"<<Log::endl;
      GriddedData griddedData(Ellipsoid(a,f), points, areas, field);
      writeFileGriddedData(fileNameGrid, griddedData);
      MiscGriddedData::printStatistics(griddedData);
    } // if(Parallel::isMaster(comm))
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
