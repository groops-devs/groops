/***********************************************/
/**
* @file griddedDatatimeSeries2GriddedData.cpp
*
* @brief Write a griddedDataTimeSeries as griddedData for each epoch.
*
* @author Torsten Mayer-Guerr
* @date 2021-09-11
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Read a \configFile{inputfileGriddedDataTimeSeries}{griddedDataTimeSeries}
and write for each epoch a \file{gridded data file}{griddedData} where
the \config{variableLoopTime} and \config{variableLoopIndex} are expanded for
each point of the given \configClass{timeSeries}{timeSeriesType}
to create the file name for this epoch (see \reference{text parser}{general.parser:text}).

If \configClass{timeSeries}{timeSeriesType} is not set
the temporal nodal points from the inputfile are used.

See also \program{GriddedData2GriddedDataTimeSeries}.
)";


/***********************************************/

#include "programs/program.h"
#include "files/fileGriddedData.h"
#include "files/fileGriddedDataTimeSeries.h"
#include "classes/timeSeries/timeSeries.h"
#include "misc/miscGriddedData.h"

/***** CLASS ***********************************/

/** @brief Write a griddedDataTimeSeries as griddedData for each epoch.
* @ingroup programsGroup */
class GriddedDataTimeSeries2GriddedData
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GriddedDataTimeSeries2GriddedData, SINGLEPROCESS, "Write a griddedDataTimeSeries as griddedData for each epoch", Grid, TimeSeries)

/***********************************************/

void GriddedDataTimeSeries2GriddedData::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName      fileNameOut, fileNameIn;
    std::string   nameTime, nameIndex, nameCount;
    TimeSeriesPtr timeSeries;

    readConfig(config, "outputfilesGriddedData",         fileNameOut, Config::MUSTSET,  "grid_{loopTime:%y-%m}.dat", "for each epoch");
    readConfig(config, "variableLoopTime",               nameTime,    Config::OPTIONAL, "loopTime", "variable with time of each epoch");
    readConfig(config, "variableLoopIndex",              nameIndex,   Config::OPTIONAL, "",         "variable with index of current epoch (starts with zero)");
    readConfig(config, "variableLoopCount",              nameCount,   Config::OPTIONAL, "",         "variable with total number of epochs");
    readConfig(config, "inputfileGriddedDataTimeSeries", fileNameIn,  Config::MUSTSET,  "",         "");
    readConfig(config, "timeSeries",                     timeSeries,  Config::OPTIONAL, "",         "otherwise times from inputfile are used");
    if(isCreateSchema(config)) return;

    logStatus<<"read gridded data time series <"<<fileNameIn<<">"<<Log::endl;
    InFileGriddedDataTimeSeries file(fileNameIn);
    GriddedData grid = file.grid();
    MiscGriddedData::printStatistics(grid);
    grid.values = std::vector<std::vector<Double>>(file.dataCount(), std::vector<Double>(grid.points.size()));
    std::vector<Time> times = file.times();
    if(timeSeries)
      times = timeSeries->times();

    VariableList varList;
    if(!nameTime.empty())  varList.undefineVariable(nameTime);
    if(!nameIndex.empty()) varList.undefineVariable(nameIndex);
    if(!nameCount.empty()) varList.setVariable(nameCount, times.size());
    for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
    {
      if(!nameTime.empty())  varList.setVariable(nameTime, times.at(idEpoch).mjd());
      if(!nameIndex.empty()) varList.setVariable(nameIndex, idEpoch);
      logStatus<<"write gridded data <"<<fileNameOut(varList)<<">"<<Log::endl;
      Matrix data = file.data(times.at(idEpoch));
      for(UInt i=0; i<grid.points.size(); i++)
        for(UInt k=0; k<data.columns(); k++)
          grid.values.at(k).at(i) = data(i, k);
      writeFileGriddedData(fileNameOut(varList), grid);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
