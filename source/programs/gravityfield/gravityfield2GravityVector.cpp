/***********************************************/
/**
* @file gravityfield2GravityVector.cpp
*
* @brief Gravity vector on a grid.
*
* @author Torsten Mayer-Guerr
* @date 2025-03-29
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program computes gravity vectors from \configClass{gravityfield}{gravityfieldType}
on a \configClass{grid}{gridType} in a global terrestrial reference frame
or alternatively in a local elliposidal frame (north, east, up) if \config{localReferenceFrame} is set.
In \configFile{outputfileGriddedData}{griddedData} the values $[gx, gy, gz]$
will be saved together with points expressed as ellipsoidal coordinates
(longitude, latitude, height) based on a reference ellipsoid with parameters \config{R} and \config{inverseFlattening}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileGriddedData.h"
#include "classes/grid/grid.h"
#include "classes/gravityfield/gravityfield.h"
#include "misc/miscGriddedData.h"

/***** CLASS ***********************************/

/** @brief Gravity vector on a grid.
* @ingroup programsGroup */
class Gravityfield2GravityVector
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Gravityfield2GravityVector, PARALLEL, "Gravity vector from gravity field", Gravityfield)

/***********************************************/

void Gravityfield2GravityVector::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName        fileNameGrid;
    GridPtr         grid;
    GravityfieldPtr gravityfield;
    Time            time;
    Bool            useLocalFrame;
    Double          a, f;

    readConfig(config, "outputfileGriddedData", fileNameGrid,   Config::MUSTSET,  "", "gx, gy, gz");
    readConfig(config, "grid",                  grid,           Config::MUSTSET,  "", "");
    readConfig(config, "gravityfield",          gravityfield,   Config::MUSTSET,  "", "");
    readConfig(config, "localReferenceFrame",   useLocalFrame,  Config::OPTIONAL, "0", "local left handed reference frame (north, east, up)");
    readConfig(config, "time",                  time,           Config::OPTIONAL, "", "at this time the gravity field will be evaluated");
    readConfig(config, "R",                     a,              Config::DEFAULT,  STRING_DEFAULT_GRS80_a, "reference radius for ellipsoidal coordinates on output");
    readConfig(config, "inverseFlattening",     f,              Config::DEFAULT,  STRING_DEFAULT_GRS80_f, "reference flattening for ellipsoidal coordinates on output, 0: spherical coordinates");
    if(isCreateSchema(config)) return;

    // create grid
    // -----------
    std::vector<Vector3d> points = grid->points();
    std::vector<Double>   areas  = grid->areas();
    Ellipsoid ellipsoid(a,f);

    // Compute gravity
    // ---------------
    logStatus<<"compute gravity"<<Log::endl;
    std::vector<Vector3d> gravity(points.size());
    Parallel::forEach(gravity, [&](UInt i){return gravityfield->gravity(time, points.at(i));}, comm);

    if(Parallel::isMaster(comm))
    {
      std::vector<std::vector<Double>> field(3, std::vector<Double>(points.size()));
      for(UInt i=0; i<points.size(); i++)
      {
        if(useLocalFrame)
          gravity.at(i) = localNorthEastUp(points.at(i), ellipsoid).inverseTransform(gravity.at(i));
        field.at(0).at(i) = gravity.at(i).x();
        field.at(1).at(i) = gravity.at(i).y();
        field.at(2).at(i) = gravity.at(i).z();
      }

      // write results
      // -------------
      logStatus<<"save values <"<<fileNameGrid<<">"<<Log::endl;
      GriddedData griddedData(ellipsoid, points, areas, field);
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
