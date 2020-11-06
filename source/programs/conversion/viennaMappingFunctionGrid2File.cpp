/***********************************************/
/**
* @file viennaMappingFunctionGrid2File.cpp
*
* @brief Converts VMF grid time series to one file.
*
* @author Torsten Mayer-Guerr
* @date 2010-04-03
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program converts the gridded time series of the Vienna Mapping Functions (VMF) into
the \file{GROOPS file format}{griddedDataTimeSeries}.

Gridded VMF data is available at: \url{https://vmf.geo.tuwien.ac.at/trop_products/GRID/}
)";

/***********************************************/

#include "programs/program.h"
#include "base/griddedData.h"
#include "inputOutput/file.h"
#include "files/fileGriddedDataTimeSeries.h"
#include "classes/timeSeries/timeSeries.h"

/***** CLASS ***********************************/

/** @brief Converts VMF grid time series to one file.
* @ingroup programsConversionGroup */
class ViennaMappingFunctionGrid2File
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(ViennaMappingFunctionGrid2File, SINGLEPROCESS, "converts VMF grid time series to one file.", Conversion)

/***********************************************/

void ViennaMappingFunctionGrid2File::run(Config &config)
{
  try
  {
    FileName              fileNameCoefficients;
    std::vector<FileName> fileNamesIn;
    TimeSeriesPtr         timeSeries;
    Angle                 deltaL, deltaB;
    Bool                  isCellRegistered;

    readConfig(config, "outputfileVmfCoefficients", fileNameCoefficients, Config::OPTIONAL, "",  "");
    readConfig(config, "inputfile",                 fileNamesIn,          Config::MUSTSET,  "",  "files must be given for each point in time");
    readConfig(config, "timeSeries",                timeSeries,           Config::MUSTSET,  "",  "times of input files");
    readConfig(config, "deltaLambda",               deltaL,               Config::DEFAULT,  "1", "[deg] sampling in longitude");
    readConfig(config, "deltaPhi",                  deltaB,               Config::DEFAULT,  "1", "[deg] sampling in latitude");
    readConfig(config, "isCellRegistered",          isCellRegistered,     Config::DEFAULT,  "1", "grid points represent cells (VMF3), not grid corners (VMF1)");
    if(isCreateSchema(config)) return;

    // ======================================================

    const std::vector<Time> times = timeSeries->times();
    if(times.size()!=fileNamesIn.size())
      throw(Exception("fileCount("+fileNamesIn.size()%"%i) != timeCount("s+times.size()%"%i)"s));

    std::vector<Angle> longitudes(static_cast<UInt>(std::round(2*PI/deltaL)));                           //!< Longitude (columns)
    for(UInt i=0; i<longitudes.size(); i++)
      longitudes.at(i) = Angle((i+(isCellRegistered ? 0.5 : 0))*deltaL);

    std::vector<Angle> latitudes(static_cast<UInt>(std::round(PI/deltaB))+(isCellRegistered ? 0 : 1));  //!< Latitude (rows)
    for(UInt i=0; i<latitudes.size(); i++)
      latitudes.at(i) = Angle(PI/2-(i+(isCellRegistered ? 0.5 : 0))*deltaB);

    GriddedDataRectangular gridRectangular;
    gridRectangular.longitudes = longitudes;
    gridRectangular.latitudes  = latitudes;
    gridRectangular.heights.resize(latitudes.size(), 0.);
    GriddedData grid;
    gridRectangular.convert(grid);

    // ======================================================

    logStatus<<"read coefficients from files"<<Log::endl;
    std::vector<Matrix> data(times.size(), Matrix(grid.points.size(), 8));
    logTimerStart;
    for(UInt i=0; i<fileNamesIn.size(); i++)
    {
      logTimerLoop(i, fileNamesIn.size());
      InFile file(fileNamesIn.at(i));

      std::string line;
      Double lon, lat;
      for(UInt k=0; k<grid.points.size(); k++)
      {
        while(std::getline(file, line) && line.size() && (line.at(0) == '!' || line.at(0) == '#')); // skip header/comment lines
        std::stringstream ss(line);
        ss>>lat>>lon;
        ss>>data.at(i)(k, 0)>>data.at(i)(k, 1)>>data.at(i)(k, 2)>>data.at(i)(k, 3); // ah, aw, zhd, zwd
        ss>>data.at(i)(k, 4)>>data.at(i)(k, 5)>>data.at(i)(k, 6)>>data.at(i)(k, 7); // gradients
      }
      data.at(i).column(4, 4) *= 0.001; // gradients mm --> m
    }
    logTimerLoopEnd(fileNamesIn.size());

    // ======================================================

    if(!fileNameCoefficients.empty())
    {
      logStatus<<"write VMF coefficients to <"<<fileNameCoefficients<<">"<<Log::endl;
      writeFileGriddedDataTimeSeries(fileNameCoefficients, 1, times, grid, data);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/
