/***********************************************/
/**
* @file griddedData2TimeSeries.cpp
*
* @brief Write time series of gridded data as time series file.
*
* @author Torsten Mayer-Guerr
* @date 2019-01-29
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Write a series of \configFile{inputfileGriddedData}{griddedData} with the corresponding
\configClass{timeSeries}{timeSeriesType} as a single time series file
(\file{instrument}{instrument}, MISCVALUES).

If \config{groupDataByPoints} is true the \config{outputfileTimeSeries} starts
for each epoch with all data (\verb|data0|, \verb|data1|\ldots) for the first point,
followed by all data of the second point and so on.
If \config{groupDataByPoints} is false, the file starts with \verb|data0|
for all points, followed by all \verb|data1| and so on.

This enables the use of all instrument programs like \program{InstrumentFilter} or
\program{InstrumentDetrend} to analyze and manipulate time series of gridded data.

See also \program{TimeSeries2GriddedData}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "files/fileGriddedData.h"
#include "classes/timeSeries/timeSeries.h"
#include "misc/miscGriddedData.h"

/***** CLASS ***********************************/

/** @brief Write time series of gridded data as time series file.
* @ingroup programsGroup */
class GriddedData2TimeSeries
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GriddedData2TimeSeries, SINGLEPROCESS, "Write time series of gridded data as time series file", Grid, TimeSeries)

/***********************************************/

void GriddedData2TimeSeries::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName              fileNameInstrument;
    std::vector<FileName> fileNamesGrid;
    TimeSeriesPtr         timeSeries;
    Bool                  groupData;

    readConfig(config, "outputfileTimeSeries", fileNameInstrument, Config::MUSTSET, "",  "each epoch: multiple data for points (MISCVALUES)");
    readConfig(config, "inputfileGriddedData", fileNamesGrid,      Config::MUSTSET, "",  "file count must agree with number of times");
    readConfig(config, "timeSeries",           timeSeries,         Config::MUSTSET, "",  "");
    readConfig(config, "groupDataByPoints",    groupData,          Config::DEFAULT, "1", "multiple data are given point by point, otherwise: data0 for all points, followed by all data1");
    if(isCreateSchema(config)) return;

    const std::vector<Time> times = timeSeries->times();
    if(times.size() != fileNamesGrid.size())
      throw(Exception("number of files ("+fileNamesGrid.size()%"%i) must agree with number of times ("s+times.size()%"%i)"s));

    MiscValuesArc arc;
    for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
    {
      try
      {
        logStatus<<"read gridded data <"<<fileNamesGrid.at(idEpoch)<<">"<<Log::endl;
        GriddedData grid;
        readFileGriddedData(fileNamesGrid.at(idEpoch), grid);

        MiscValuesEpoch epoch(grid.points.size() * grid.values.size());
        epoch.time = times.at(idEpoch);
        for(UInt k=0; k<grid.values.size(); k++)
          for(UInt i=0; i<grid.values.at(k).size(); i++)
            epoch.values(groupData ? (i*grid.values.size()+k) : (i+grid.points.size()*k)) = grid.values.at(k).at(i);
        arc.push_back(epoch);
      }
      catch(std::exception &e)
      {
        logWarning<<e.what()<<Log::endl;
      }
    }

    logStatus<<"write time series to file <"<<fileNameInstrument<<">"<<Log::endl;
    InstrumentFile::write(fileNameInstrument, arc);
    Arc::printStatistics(arc);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
