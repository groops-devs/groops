/***********************************************/
/**
* @file gravityfieldVariancesPropagation2GriddedData.cpp
*
* @brief Standard deviations of values of a gravity field on a grid.
*
* @author Torsten Mayer-Guerr
* @date 2008-11-06
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program propagates variance-covariance matrix of a \configClass{gravityfield}{gravityfieldType}
evaluated at \config{time} to the points of a \configClass{grid}{gridType} in terms of the functional
given by \configClass{kernel}{kernelType}.
The resulting \file{outputfileGriddedData}{griddedData} contains the standard deviations of the grid
points.

See also \program{Gravityfield2GridCovarianceMatrix}, \program{GravityfieldCovariancesPropagation2GriddedData}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileGriddedData.h"
#include "classes/grid/grid.h"
#include "classes/kernel/kernel.h"
#include "classes/gravityfield/gravityfield.h"
#include "misc/miscGriddedData.h"

/***** CLASS ***********************************/

/** @brief Standard deviations of values of a gravity field on a grid.
* @ingroup programsGroup */
class GravityfieldVariancesPropagation2GriddedData
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GravityfieldVariancesPropagation2GriddedData, PARALLEL, "standard deviations of values of a gravity field on a grid", Gravityfield, Grid, Covariance)

/***********************************************/

void GravityfieldVariancesPropagation2GriddedData::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName        fileNameGrid;
    GridPtr         grid;
    KernelPtr       kernel;
    GravityfieldPtr gravityfield;
    Time            time;
    Double          a, f;

    readConfig(config, "outputfileGriddedData", fileNameGrid, Config::MUSTSET,  "", "standard deviation at each grid point");
    readConfig(config, "grid",                  grid,         Config::MUSTSET,  "", "");
    readConfig(config, "kernel",                kernel,       Config::MUSTSET,  "", "functional");
    readConfig(config, "gravityfield",          gravityfield, Config::MUSTSET,  "", "");
    readConfig(config, "time",                  time,         Config::OPTIONAL, "", "at this time the gravity field will be evaluated");
    readConfig(config, "R",                     a,            Config::DEFAULT,  STRING_DEFAULT_GRS80_a, "reference radius for ellipsoidal coordinates on output");
    readConfig(config, "inverseFlattening",     f,            Config::DEFAULT,  STRING_DEFAULT_GRS80_f, "reference flattening for ellipsoidal coordinates on output, 0: spherical coordinates");
    if(isCreateSchema(config)) return;

    // Compute standard deviations
    // ---------------------------
    logStatus<<"calculate standard deviations on grid"<<Log::endl;
    std::vector<Vector3d> points = grid->points();
    std::vector<Double>   areas  = grid->areas();
    std::vector<Double>   field(points.size());
    Parallel::forEach(field, [&](UInt i) {return sqrt(gravityfield->variance(time, points.at(i), *kernel));}, comm);

    if(Parallel::isMaster(comm))
    {
      logStatus<<"save values to file <"<<fileNameGrid<<">"<<Log::endl;
      GriddedData griddedData(Ellipsoid(a,f), points, areas, {field});
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
