/***********************************************/
/**
* @file gravityfield2GriddedData.cpp
*
* @brief Values of a gravity field on a grid.
*
* @author Torsten Mayer-Guerr
* @date 2002-03-19
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program computes values of a \configClass{gravityfield}{gravityfieldType} on a given \configClass{grid}{gridType}.
The type of value (e.g gravity anomalies or geoid heights) can be choosen with \configClass{kernel}{kernelType}.
If a time is given the gravity field will be evaluated at this point of time otherwise only the static part will be used.
The values will be saved together with points expressed as ellipsoidal coordinates (longitude, latitude, height)
based on a reference ellipsoid with parameters \config{R} and \config{inverseFlattening}.
To speed up the computation the gravity field can be converted to spherical harmonics before the computation
with \config{convertToHarmonics}.

To visualize the results use \program{PlotMap}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileGriddedData.h"
#include "classes/grid/grid.h"
#include "classes/kernel/kernel.h"
#include "classes/gravityfield/gravityfield.h"
#include "misc/miscGriddedData.h"

/***** CLASS ***********************************/

/** @brief Values of a gravity field on a grid.
* @ingroup programsGroup */
class Gravityfield2GriddedData
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Gravityfield2GriddedData, PARALLEL, "values of a gravity field on a grid.", Gravityfield, Grid, Covariance)

/***********************************************/

void Gravityfield2GriddedData::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName          fileNameGrid;
    GridPtr           grid;
    KernelPtr         kernel;
    GravityfieldPtr   gravityfield;
    Time              time;
    Double            a, f;
    Bool              convertToHarmonics;

    readConfig(config, "outputfileGriddedData", fileNameGrid,       Config::MUSTSET,  "", "");
    readConfig(config, "grid",                  grid,               Config::MUSTSET,  "", "");
    readConfig(config, "kernel",                kernel,             Config::MUSTSET,  "", "");
    readConfig(config, "gravityfield",          gravityfield,       Config::MUSTSET,  "", "");
    readConfig(config, "convertToHarmonics",    convertToHarmonics, Config::DEFAULT,   "1", "gravityfield is converted to spherical harmonics before evaluation, may accelerate the computation");
    readConfig(config, "time",                  time,               Config::OPTIONAL,  "", "at this time the gravity field will be evaluated");
    readConfig(config, "R",                     a,                  Config::DEFAULT,   STRING_DEFAULT_GRS80_a, "reference radius for ellipsoidal coordinates on output");
    readConfig(config, "inverseFlattening",     f,                  Config::DEFAULT,   STRING_DEFAULT_GRS80_f, "reference flattening for ellipsoidal coordinates on output, 0: spherical coordinates");
    if(isCreateSchema(config)) return;

    // Create values on grid
    // ---------------------
    logStatus<<"create values on grid"<<Log::endl;
    std::vector<Double> field(grid->points().size());
    if(!convertToHarmonics) // All representations, all point distributions
      Parallel::forEach(field, [&](UInt i){return gravityfield->field(time, grid->points().at(i), *kernel);}, comm);
    else // fast version
      field = MiscGriddedData::synthesisSphericalHarmonics(gravityfield->sphericalHarmonics(time), grid->points(), kernel, comm);

    if(Parallel::isMaster(comm))
    {
      logStatus<<"save values to file <"<<fileNameGrid<<">"<<Log::endl;
      GriddedData griddedData(Ellipsoid(a,f), grid->points(), grid->areas(), {field});
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
