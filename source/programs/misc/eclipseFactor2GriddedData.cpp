/***********************************************/
/**
* @file eclipseFactor2GriddedData.cpp
*
* @brief Eclipse factor on a grid.
*
* @author Andreas Strasser
* @date 2023-10-13
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program converts the output of a \configClass{eclipse}{eclipseType} model on a given
\configClass{grid}{gridType}. The time for the evaluation can be specified in \config{time}.
The values will be saved together with points expressed as ellipsoidal coordinates
(longitude, latitude, height) based on a reference ellipsoid with parameters \config{R}
and \config{inverseFlattening}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileGriddedData.h"
#include "classes/grid/grid.h"
#include "classes/eclipse/eclipse.h"
#include "classes/ephemerides/ephemerides.h"
#include "misc/miscGriddedData.h"
#include "classes/earthRotation/earthRotation.h"

/***** CLASS ***********************************/

/** @brief Eclipse factor on a grid.
* @ingroup programsGroup */
class EclipseFactor2GriddedData
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(EclipseFactor2GriddedData, PARALLEL, "eclipse factor on a grid.", Misc, Grid)

/***********************************************/

void EclipseFactor2GriddedData::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName         fileNameGrid;
    GridPtr          grid;
    EclipsePtr       eclipse;
    EphemeridesPtr   ephemerides;
    EarthRotationPtr earthRotation;
    Time             time;
    Double           a, f;

    readConfig(config, "outputfileGriddedData", fileNameGrid,  Config::MUSTSET,  "",  "eclipse factor");
    readConfig(config, "grid",                  grid,          Config::MUSTSET,  "",  "");
    readConfig(config, "eclipse",               eclipse,       Config::MUSTSET,  "",  "");
    readConfig(config, "ephemerides",           ephemerides,   Config::MUSTSET,  "",  "");
    readConfig(config, "earthRotation",         earthRotation, Config::OPTIONAL, "file", "");
    readConfig(config, "time",                  time,          Config::MUSTSET,  "",  "");
    readConfig(config, "R",                     a,             Config::DEFAULT,  STRING_DEFAULT_GRS80_a, "reference radius for ellipsoidal coordinates on output");
    readConfig(config, "inverseFlattening",     f,             Config::DEFAULT,  STRING_DEFAULT_GRS80_f, "reference flattening for ellipsoidal coordinates on output, 0: spherical coordinates");
    if(isCreateSchema(config)) return;

    // compute
    // -------
    logStatus<<"computing eclipse factor"<<Log::endl;
    Rotary3d rotation;
    if(earthRotation)
      rotation = earthRotation->rotaryMatrix(time);
    std::vector<Vector3d> points = grid->points();
    std::vector<Double>   values(points.size());

    Parallel::forEach(values, [&](UInt i)
    {
      return eclipse->factor(time, rotation.inverseRotate(points.at(i)), ephemerides);
    }, comm);

    if(Parallel::isMaster(comm))
    {
      logStatus<<"save values to file <"<<fileNameGrid<<">"<<Log::endl;
      GriddedData griddedData(Ellipsoid(a, f), points, grid->areas(), {values});
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
