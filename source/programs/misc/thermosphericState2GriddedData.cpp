/***********************************************/
/**
* @file thermosphericState2GriddedData.cpp
*
* @brief Thermospheric state values.
**
* @author Torsten Mayer-Guerr
* @date 2020-07-05
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program converts the output (neutral mass density,temperature) of an empirical thermosphere model (e.g. JB2008) on a given \configClass{grid}{gridType}. 
Additionally, also the thermospheric winds estimated by using the horizontal wind model HWM 2014 can be assessed.  
The time for the evaluation can be specified in \config{time}. The values will be saved together with points expressed as ellipsoidal coordinates
(longitude, latitude, height) based on a reference ellipsoid with parameters \config{R} and \config{inverseFlattening}.    
\fig{!hb}{1.0}{thermosphericState2GriddedData}{fig:thermosphericState2GriddedData}{JB2008 model in 300 km height at 2003-07-01 12:00.}
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileGriddedData.h"
#include "classes/grid/grid.h"
#include "classes/thermosphere/thermosphere.h"
#include "misc/miscGriddedData.h"

/***** CLASS ***********************************/

/** @brief Thermospheric state values.
* @ingroup programsGroup */
class ThermosphericState2GriddedData
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(ThermosphericState2GriddedData, PARALLEL, "thermospheric state values", Misc, Grid)

/***********************************************/

void ThermosphericState2GriddedData::run(Config &config)
{
  try
  {
    FileName        fileNameGrid;
    ThermospherePtr thermosphere;
    GridPtr         grid;
    Time            time;
    Bool            useLocalFrame;
    Double          a, f;

    readConfig(config, "outputfileGriddedData", fileNameGrid,  Config::MUSTSET,  "", "density [kg/m**3], temperature [K], wind (x, y, z) [m/s**2]");
    readConfig(config, "thermosphere",          thermosphere,  Config::MUSTSET,  "", "");
    readConfig(config, "grid",                  grid,          Config::MUSTSET,  "", "");
    readConfig(config, "time",                  time,          Config::MUSTSET,  "", "");
    readConfig(config, "localReferenceFrame",   useLocalFrame, Config::OPTIONAL, "1", "wind in local north, east, up, otherwise global terrestrial");
    readConfig(config, "R",                     a,             Config::DEFAULT,  STRING_DEFAULT_GRS80_a, "reference radius for ellipsoidal coordinates on output");
    readConfig(config, "inverseFlattening",     f,             Config::DEFAULT,  STRING_DEFAULT_GRS80_f, "reference flattening for ellipsoidal coordinates on output, 0: spherical coordinates");
    if(isCreateSchema(config)) return;

    // compute
    // -------
    logStatus<<"compute thermosphere"<<Log::endl;
    Ellipsoid             ellipsoid(a, f);
    std::vector<Vector3d> points = grid->points();
    std::vector<Double>   areas  = grid->areas();
    std::vector<Vector>   values(points.size());
    Parallel::forEach(values, [&](UInt i)
    {
      Double   density, temperature;
      Vector3d wind;
      thermosphere->state(time, points.at(i), density, temperature, wind);
      if(useLocalFrame)
        wind = localNorthEastUp(points.at(i), ellipsoid).inverseTransform(wind);
      return Vector({density, temperature, wind.x(), wind.y(), wind.z()});
    });

    if(Parallel::isMaster())
    {
      // convert
      std::vector<std::vector<Double>> field(5, std::vector<Double>(points.size()));
      for(UInt i=0; i<points.size(); i++)
        for(UInt k=0; k<field.size(); k++)
          field.at(k).at(i) = values.at(i).at(k);

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
