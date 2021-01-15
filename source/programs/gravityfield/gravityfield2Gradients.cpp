/***********************************************/
/**
* @file gravityfield2Gradients.cpp
*
* @brief Gradients on a grid.
*
* @author Daniel Rieser
* @date 2015-05-11
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program computes gravity gradients from \configClass{gravityfield}{gravityfieldType}
on a \configClass{grid}{gridType} in a global terrestrial reference frame
or alternatively in a local elliposidal frame (north, east, up) if \config{localReferenceFrame} is set.
In \configFile{outputfileGriddedData}{griddedData} the values $[Vxx, Vyy, Vzz, Vxy, Vxz, Vyz]$
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

/** @brief Deflections of the vertical on a grid.
* @ingroup programsGroup */
class Gravityfield2Gradients
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Gravityfield2Gradients, PARALLEL, "Gradients from gravity field", Gravityfield)

/***********************************************/

void Gravityfield2Gradients::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName        fileNameGrid;
    GridPtr         grid;
    GravityfieldPtr gravityfield;
    Time            time;
    Bool            useLocalFrame;
    Double          a, f;

    readConfig(config, "outputfileGriddedData", fileNameGrid,   Config::MUSTSET,  "", "Vxx Vyy Vzz Vxy Vxz Vyz");
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

    // Compute gradients
    // ------------------
    logStatus<<"compute gradients"<<Log::endl;
    std::vector<Tensor3d> gradients(points.size());
    Parallel::forEach(gradients, [&](UInt i){return gravityfield->gravityGradient(time, points.at(i));}, comm);

    if(Parallel::isMaster(comm))
    {
      // convert
      // -------
      std::vector<std::vector<Double>> field(6, std::vector<Double>(points.size()));
      for(UInt i=0; i<points.size(); i++)
      {
        if(useLocalFrame)
          gradients.at(i) = localNorthEastUp(points.at(i), ellipsoid).inverseTransform(gradients.at(i));
        field.at(0).at(i) = gradients.at(i).xx();
        field.at(1).at(i) = gradients.at(i).yy();
        field.at(2).at(i) = gradients.at(i).zz();
        field.at(3).at(i) = gradients.at(i).xy();
        field.at(4).at(i) = gradients.at(i).xz();
        field.at(5).at(i) = gradients.at(i).yz();
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
