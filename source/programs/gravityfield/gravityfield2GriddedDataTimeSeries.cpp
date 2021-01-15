/***********************************************/
/**
* @file gravityfield2GriddedDataTimeSeries.cpp
*
* @brief Time series of gridded gravity fields.
*
* @author Torsten Mayer-Guerr
* @date 2020-07-20
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program computes values of a \configClass{gravityfield}{gravityfieldType} on a given \configClass{grid}{gridType}
for each time step of \configClass{timeSeries}{timeSeriesType}.
The type of value (e.g gravity anomalies or geoid heights) can be choosen with \configClass{kernel}{kernelType}.
To speed up the computation the gravity field can be converted to spherical harmonics before the computation
with \config{convertToHarmonics}.
The \configFile{outputfileTimeSeries}{instrument} is an instrument (MISCVALUES) file with a data column
for each grid point per epoch.

This program enables the use of all instrument programs like \program{InstrumentFilter},
\program{InstrumentArcStatistics} or \program{InstrumentDetrend} to analyze and manipulate time series of gridded data.

See also \program{TimeSeries2GriddedData}, \program{Gravityfield2GriddedData}
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "classes/grid/grid.h"
#include "classes/kernel/kernel.h"
#include "classes/gravityfield/gravityfield.h"
#include "classes/timeSeries/timeSeries.h"
#include "misc/miscGriddedData.h"

/***** CLASS ***********************************/

/** @brief Time series of gridded gravity fields.
* @ingroup programsGroup */
class Gravityfield2GriddedDataTimeSeries
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Gravityfield2GriddedDataTimeSeries, PARALLEL, "time series of gridded gravity fields.", Gravityfield, Grid, TimeSeries)

/***********************************************/

void Gravityfield2GriddedDataTimeSeries::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName        fileNameOut;
    GridPtr         grid;
    KernelPtr       kernel;
    GravityfieldPtr gravityfield;
    TimeSeriesPtr   timeSeries;
    Bool            convertToHarmonics;

    readConfig(config, "outputfileTimeSeries", fileNameOut,        Config::MUSTSET, "",  "each epoch: data of grid points (MISCVALUES)");
    readConfig(config, "grid",                 grid,               Config::MUSTSET, "",  "");
    readConfig(config, "kernel",               kernel,             Config::MUSTSET, "",  "");
    readConfig(config, "gravityfield",         gravityfield,       Config::MUSTSET, "",  "");
    readConfig(config, "convertToHarmonics",   convertToHarmonics, Config::DEFAULT, "1", "gravityfield is converted to spherical harmonics before evaluation, may accelerate the computation");
    readConfig(config, "timeSeries",           timeSeries,         Config::MUSTSET, "",  "");
    if(isCreateSchema(config)) return;

    // create grid
    // -----------
    const std::vector<Vector3d> points = grid->points();
    const std::vector<Double>   areas  = grid->areas();
    const std::vector<Time>     times  = timeSeries->times();

    // Create values on grid
    // ---------------------
    logStatus<<"create values on grid for each time"<<Log::endl;
    Matrix A(times.size(), 1+points.size()); // one time column + data
    Parallel::forEach(times.size(), [&](UInt i)
    {
      if(!convertToHarmonics) // All representations, all point distributions
      {
        for(UInt k=0; k<points.size(); k++)
          A(i, 1+k) = gravityfield->field(times.at(i), points.at(k), *kernel);
      }
      else // fast version
      {
        std::vector<Double> field = MiscGriddedData::synthesisSphericalHarmonics(gravityfield->sphericalHarmonics(times.at(i)), points, kernel, Parallel::selfCommunicator(), FALSE/*timing*/);
        copy(Vector(field).trans(), A.slice(i, 1, 1, field.size()));
      }
    }, comm);
    Parallel::reduceSum(A, 0, comm);

    // Write results
    // -------------
    if(Parallel::isMaster(comm))
    {
      logStatus<<"write time series to file <"<<fileNameOut<<">"<<Log::endl;
      InstrumentFile::write(fileNameOut, Arc(times, A));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
