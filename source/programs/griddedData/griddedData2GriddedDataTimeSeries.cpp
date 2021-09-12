/***********************************************/
/**
* @file griddedData2GriddedDataTimeSeries.cpp
*
* @brief Write time series of gridded data as gridded data time series file.
*
* @author Torsten Mayer-Guerr
* @date 2021-09-11
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Write a series of \configFile{inputfileGriddedData}{griddedData}
with the corresponding \configClass{timeSeries}{timeSeriesType}
as a single \file{gridded data time series file}{griddedDataTimeSeries}.
The \config{splineDegree} defines the possible temporal interpolation of data in the output file.
For a file with spline degree 0 (temporal block means) the time intervals
in which the grids are valid are defined between adjacent points in time.
Therefore one more point in time is needed than the number of input grid files for degree 0.

See also \program{GriddedDataTimeSeries2GriddedData}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileGriddedData.h"
#include "files/fileGriddedDataTimeSeries.h"
#include "classes/timeSeries/timeSeries.h"

/***** CLASS ***********************************/

/** @brief Write time series of gridded data as gridded data time series file.
* @ingroup programsGroup */
class GriddedData2GriddedDataTimeSeries
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GriddedData2GriddedDataTimeSeries, SINGLEPROCESS, "Write time series of gridded data as gridded data time series file", Grid, TimeSeries)

/***********************************************/

void GriddedData2GriddedDataTimeSeries::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName              fileNameOut;
    std::vector<FileName> fileNamesGrid;
    TimeSeriesPtr         timeSeries;
    UInt                  splineDegree;

    readConfig(config, "outputfileGriddedDataTimeSeries", fileNameOut,   Config::MUSTSET, "",  "");
    readConfig(config, "inputfileGriddedData",            fileNamesGrid, Config::MUSTSET, "",  "file count must agree with number of times+splineDegre-1");
    readConfig(config, "timeSeries",                      timeSeries,    Config::MUSTSET, "",  "");
    readConfig(config, "splineDegree",                    splineDegree,  Config::DEFAULT, "1", "degree of splines");
    if(isCreateSchema(config)) return;

    const std::vector<Time> times = timeSeries->times();
    if(times.size()+splineDegree-1 != fileNamesGrid.size())
      throw(Exception("fileCount("+fileNamesGrid.size()%"%i) != timeCount("s+times.size()%"%i)-1+splineDegree"s));

    GriddedData         grid;
    std::vector<Matrix> data;
    for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
    {
      try
      {
        logStatus<<"read gridded data <"<<fileNamesGrid.at(idEpoch)<<">"<<Log::endl;
        readFileGriddedData(fileNamesGrid.at(idEpoch), grid);
        if(!data.size())
          data.resize(times.size(), Matrix(grid.points.size(), grid.values.size()));
        for(UInt k=0; k<grid.values.size(); k++)
          for(UInt i=0; i<grid.values.at(k).size(); i++)
            data.at(idEpoch)(i,k) = grid.values.at(k).at(i);
      }
      catch(std::exception &e)
      {
        logWarning<<e.what()<<Log::endl;
      }
    }

    logStatus<<"write time series to file <"<<fileNameOut<<">"<<Log::endl;
    grid.values.clear();
    writeFileGriddedDataTimeSeries(fileNameOut, splineDegree, times, grid, data);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
