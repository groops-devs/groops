/***********************************************/
/**
* @file gravityfield2AbsoluteGravity.cpp
*
* @brief Absolute gravity values on a grid.
*
* @author Torsten Mayer-Guerr
* @date 2011-10-29
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program computes the absolute value of gravity $\left\lVert{\M g}\right\rVert$
of a \configClass{gravityfield}{gravityfieldType} on a given \configClass{grid}{gridType}.
The result is multiplicated with \config{factor}.
To get the full gravity vector in a terrestrial frame add
the centrifugal part, see \configClass{gravityfield:tides:centrifugal}{tidesType:centrifugal}.

The values will be saved together with points expressed as ellipsoidal coordinates (longitude, latitude, height)
based on a reference ellipsoid with parameters \config{R} and \config{inverseFlattening}.

It is intended to compute gravity anomalies from absolute gravity observations.
To visualize the results use \program{PlotMap}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileGriddedData.h"
#include "classes/grid/grid.h"
#include "classes/gravityfield/gravityfield.h"
#include "misc/miscGriddedData.h"

/***** CLASS ***********************************/

/** @brief Absolute gravity values on a grid.
* @ingroup programsGroup */
class Gravityfield2AbsoluteGravity
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Gravityfield2AbsoluteGravity, PARALLEL, "Absolute gravity values on a grid.", Gravityfield)

/***********************************************/

void Gravityfield2AbsoluteGravity::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName          outFileNameGrid;
    GridPtr           grid;
    GravityfieldPtr   gravityfield;
    Double            factor;
    Time              time;
    Double            a, f;

    readConfig(config, "outputfileGriddedData", outFileNameGrid, Config::MUSTSET,  "", "");
    readConfig(config, "grid",                  grid,            Config::MUSTSET,  "", "");
    readConfig(config, "gravityfield",          gravityfield,    Config::MUSTSET,  "", "");
    readConfig(config, "factor",                factor,          Config::DEFAULT,   "1.0", "the result is multiplied by this factor, set -1 to subtract the field");
    readConfig(config, "time",                  time,            Config::OPTIONAL,  "", "at this time the gravity field will be evaluated");
    readConfig(config, "R",                     a,               Config::DEFAULT,   STRING_DEFAULT_GRS80_a, "reference radius for ellipsoidal coordinates on output");
    readConfig(config, "inverseFlattening",     f,               Config::DEFAULT,   STRING_DEFAULT_GRS80_f, "reference flattening for ellipsoidal coordinates on output, 0: spherical coordinates");
    if(isCreateSchema(config)) return;

    // create grid
    // -----------
    std::vector<Vector3d> points = grid->points();
    std::vector<Double>   areas  = grid->areas();

    // Create values on grid
    // ---------------------
    logStatus<<"create values on grid"<<Log::endl;
    std::vector<Double> field(points.size());
    Parallel::forEach(field, [&](UInt i){return factor*gravityfield->gravity(time, points.at(i)).r();}, comm);

    if(Parallel::isMaster(comm))
    {
      // write results
      // -------------
      logStatus<<"save values to file <"<<outFileNameGrid<<">"<<Log::endl;
      GriddedData griddedData(Ellipsoid(a,f), points, areas, {field});
      writeFileGriddedData(outFileNameGrid, griddedData);
      MiscGriddedData::printStatistics(griddedData);
    } // Master
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
