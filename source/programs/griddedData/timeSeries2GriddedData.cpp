/***********************************************/
/**
* @file timeSeries2GriddedData.cpp
*
* @brief Write time series as gridded data for each epoch.
*
* @author Torsten Mayer-Guerr
* @date 2019-01-29
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Interpret the data columns of \configFile{inputfileTimeSeries}{instrument}
as data points of a corresponding \configClass{grid}{gridType}.

For each epoch a \file{gridded data file}{griddedData} is written where
the \config{variableLoopTime} and \config{variableLoopIndex} are expanded for
each point of the given time series to create the file name for this epoch
(see \reference{text parser}{general.parser:text}).

The number of input data columns must be a multiple of the number $n$ of grid points.
If \config{isGroupedDataByPoint} is true the \configFile{inputfileTimeSeries}{instrument} starts
with all data (\verb|data0|, \verb|data1|\ldots) for the first point, followed by all data of the second point and so on.
If \config{isGroupedDataByPoint} is false, the file starts with \verb|data0| for all points, followed by all \verb|data1| and so on.

See also \program{GriddedData2TimeSeries}.
)";


/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "files/fileGriddedData.h"
#include "classes/grid/grid.h"
#include "misc/miscGriddedData.h"

/***** CLASS ***********************************/

/** @brief Write time series as gridded data for each epoch.
* @ingroup programsGroup */
class TimeSeries2GriddedData
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(TimeSeries2GriddedData, SINGLEPROCESS, "Write time series as gridded data for each epoch", Grid, TimeSeries)

/***********************************************/

void TimeSeries2GriddedData::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName    fileNameGrid, fileNameInstrument;
    std::string nameTime, nameIndex, nameCount;
    GridPtr     gridPtr;
    Bool        groupData;
    Double      a, f;

    readConfig(config, "outputfilesGriddedData", fileNameGrid,       Config::MUSTSET,  "grid_{loopTime:%y-%m}.dat", "for each epoch");
    readConfig(config, "variableLoopTime",       nameTime,           Config::OPTIONAL, "loopTime",                  "variable with time of each epoch");
    readConfig(config, "variableLoopIndex",      nameIndex,          Config::OPTIONAL, "",                          "variable with index of current epoch (starts with zero)");
    readConfig(config, "variableLoopCount",      nameCount,          Config::OPTIONAL, "",                          "variable with total number of epochs");
    readConfig(config, "inputfileTimeSeries",    fileNameInstrument, Config::MUSTSET,  "",                          "each epoch: multiple data for points (MISCVALUES)");
    readConfig(config, "grid",                   gridPtr,            Config::MUSTSET,  "",                          "corresponding grid points");
    readConfig(config, "isDataGroupedByPoint",   groupData,          Config::DEFAULT, "1",                          "multiple data are given point by point, otherwise: first data0 for all points, followed by all data1");
    readConfig(config, "R",                      a,                  Config::DEFAULT,  STRING_DEFAULT_GRS80_a,      "reference radius for ellipsoidal coordinates on output");
    readConfig(config, "inverseFlattening",      f,                  Config::DEFAULT,  STRING_DEFAULT_GRS80_f,      "reference flattening for ellipsoidal coordinates on output, 0: spherical coordinates");
    if(isCreateSchema(config)) return;

    logStatus<<"read time series <"<<fileNameInstrument<<">"<<Log::endl;
    MiscValuesArc arc = InstrumentFile::read(fileNameInstrument);
    Arc::printStatistics(arc);

    logStatus<<"create grid"<<Log::endl;
    GriddedData grid(Ellipsoid(a,f), gridPtr->points(), gridPtr->areas(), {});

    // number of data columns
    const UInt countPoints = grid.points.size();
    const UInt countData   = arc.at(0).values.rows();
    const UInt columns     = countData/countPoints;
    if((countData < countPoints) || (countData%countPoints != 0))
      throw(Exception("number of data columns ("+countData%"%i) in time series must be a multiple of grid points ("s+countPoints%"%i)"s));
    grid.values.resize(columns, std::vector<Double>(countPoints));

    const std::vector<Time> times = arc.times();
    auto varList = config.getVarList();
    if(!nameTime.empty())  addVariable(nameTime,  varList);
    if(!nameIndex.empty()) addVariable(nameIndex, varList);
    if(!nameCount.empty()) addVariable(nameCount, times.size(), varList);

    for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
    {
      if(!nameTime.empty())  varList[nameTime]->setValue(times.at(idEpoch).mjd());
      if(!nameIndex.empty()) varList[nameIndex]->setValue(idEpoch);
      logStatus<<"write gridded data <"<<fileNameGrid(varList)<<">"<<Log::endl;

      for(UInt i=0; i<grid.points.size(); i++)
        for(UInt k=0; k<columns; k++)
          grid.values.at(k).at(i) = arc.at(idEpoch).values((groupData) ? (i*columns+k) : (i+grid.points.size()*k));

      writeFileGriddedData(fileNameGrid(varList), grid);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
