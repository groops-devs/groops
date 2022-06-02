/***********************************************/
/**
* @file magneticField2GriddedData.cpp
*
* @brief Magentic field vector.
**
* @author Torsten Mayer-Guerr
* @date 2019-05-25
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Computes x, y, z of the magentic field vector.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileGriddedData.h"
#include "classes/grid/grid.h"
#include "classes/magnetosphere/magnetosphere.h"
#include "misc/miscGriddedData.h"

/***** CLASS ***********************************/

/** @brief Magentic field vector.
* @ingroup programsGroup */
class MagneticField2GriddedData
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(MagneticField2GriddedData, PARALLEL, "magentic field vector", Misc, Grid)

/***********************************************/

void MagneticField2GriddedData::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName         fileNameGrid;
    MagnetospherePtr magnetosphere;
    GridPtr          grid;
    Time             time;
    Bool             useLocalFrame;
    Double           a, f;

    readConfig(config, "outputfileGriddedData", fileNameGrid,  Config::MUSTSET,  "", "x, y, z [Tesla = kg/A/s**2]");
    readConfig(config, "magnetosphere",         magnetosphere, Config::MUSTSET,  "", "");
    readConfig(config, "grid",                  grid,          Config::MUSTSET,  "", "");
    readConfig(config, "time",                  time,          Config::OPTIONAL, STRING_J2000, "");
    readConfig(config, "localReferenceFrame",   useLocalFrame, Config::OPTIONAL, "0", "local left handed reference frame (north, east, up)");
    readConfig(config, "R",                     a,             Config::DEFAULT,  STRING_DEFAULT_GRS80_a, "reference radius for ellipsoidal coordinates on output");
    readConfig(config, "inverseFlattening",     f,             Config::DEFAULT,  STRING_DEFAULT_GRS80_f, "reference flattening for ellipsoidal coordinates on output, 0: spherical coordinates");
    if(isCreateSchema(config)) return;

    // compute
    // -------
    logStatus<<"compute magnetic field"<<Log::endl;
    Ellipsoid             ellipsoid(a, f);
    std::vector<Vector3d> points = grid->points();
    std::vector<Double>   areas  = grid->areas();
    std::vector<Vector3d> values(points.size());
    Parallel::forEach(values, [&](UInt i) {return magnetosphere->magenticFieldVector(time, points.at(i));}, comm);

    if(Parallel::isMaster(comm))
    {
      // convert
      std::vector<std::vector<Double>> field(3, std::vector<Double>(points.size()));
      for(UInt i=0; i<points.size(); i++)
      {
        if(useLocalFrame)
          values.at(i) = localNorthEastUp(points.at(i), ellipsoid).inverseTransform(values.at(i));
        field.at(0).at(i) = values.at(i).x();
        field.at(1).at(i) = values.at(i).y();
        field.at(2).at(i) = values.at(i).z();
      }

      // write
      logStatus<<"save values to file <"<<fileNameGrid<<">"<<Log::endl;
      GriddedData griddedData(ellipsoid, points, areas, field);
      writeFileGriddedData(fileNameGrid, griddedData);
      MiscGriddedData::printStatistics(griddedData);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
