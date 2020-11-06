/***********************************************/
/**
* @file griddedData2AreaMeanTimeSeries.cpp
*
* @brief Generates a time series as mean values over an area from a list of grid files
*
* @author Andreas Kvas
* @date 2017-04-12
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program computes a time series of area mean values
in a basin represented by \configClass{border}{borderType} from a sequence of grid files.
If a file is not found, the epoch is skipped. Per default
the weighted average of all points in the given border is computed where the points are weighted by their area element.

If \config{computeMean} is set, the time average of each grid points is subtracted before the computation.
If \config{multiplyWithArea} is set, the weighted average is multiplied with the total basin area.
This is useful for computing the total mass in the basin.

The \configFile{outputfileTimeSeries}{instrument} is an instrument file, where the first columns are the
mean value each data column in the grid files, followed by the the weighted RMS
for each data column in the grid files if \config{computeRms} is set.
If the number of data columns differs between the grid files, the remaing columns are padded with zeros.

See also \program{Gravityfield2AreaMeanTimeSeries}.
)";


/***********************************************/

#include "programs/program.h"
#include "parser/dataVariables.h"
#include "files/fileGriddedData.h"
#include "files/fileInstrument.h"
#include "classes/timeSeries/timeSeries.h"
#include "classes/border/border.h"


/***** CLASS ***********************************/

/** @brief Generates a time series as mean values over an area from a list of grid files.
* @ingroup programsGroup */
class GriddedData2AreaMeanTimeSeries
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(GriddedData2AreaMeanTimeSeries, SINGLEPROCESS, "generates a time series as mean values over an area from a list of grid files", Grid, TimeSeries)

/***********************************************/

void GriddedData2AreaMeanTimeSeries::run(Config &config)
{
  try
  {
    FileName        outputName;
    std::vector<FileName> inputName;
    BorderPtr       border;
    TimeSeriesPtr   timeSeries;
    Bool            computeRms, removeMean, multiplyWithArea;

    readConfig(config, "outputfileTimeSeries", outputName,         Config::MUSTSET,  "", "");
    readConfig(config, "inputfileGriddedData", inputName ,         Config::MUSTSET, "", "");
    readConfig(config, "border",               border,             Config::MUSTSET,  "", "");
    readConfig(config, "timeSeries",           timeSeries,         Config::MUSTSET,  "", "");
    readConfig(config, "multiplyWithArea",     multiplyWithArea,   Config::DEFAULT,  "0", "multiply time series with total area (useful for mass estimates)");
    readConfig(config, "removeMean",           removeMean,         Config::DEFAULT,  "0", "remove the temporal mean of the series");
    readConfig(config, "computeRms",           computeRms,         Config::DEFAULT,  "0", "addtional rms each time step");
    if(isCreateSchema(config)) return;

    const UInt fileCount = inputName.size();
    std::vector<Time> times  = timeSeries->times();
    if(times.size() != fileCount)
      throw(Exception("Number of input files and epochs do not match (" + fileCount%"%d"s + " vs. " + times.size()%"%d"s + ")."));

    std::vector<Vector> values;
    std::vector<Vector> rms;
    std::vector<Time> gridTimes;
    UInt maxDataFields = 0;

    logTimerStart;
    for(UInt k = 0; k<times.size(); k++)
    {
      logTimerLoop(k, times.size());
      GriddedData grid;
      try
      {
        readFileGriddedData(inputName.at(k), grid);
      }
      catch(std::exception &e)
      {
        logWarning<<e.what()<<Log::endl;
        continue;
      }

      maxDataFields = std::max(maxDataFields, grid.values.size());
      values.push_back(Vector(grid.values.size()));
      rms.push_back(Vector(grid.values.size()));
      gridTimes.push_back(times.at(k));

      // find all points in border
      std::vector<UInt> pointIndex;
      for(UInt i = 0; i<grid.points.size(); i++)
        if(border->isInnerPoint(grid.points.at(i), grid.ellipsoid))
          pointIndex.push_back(i);

      // compute total area
      Double totalArea = 0.0;
      for( UInt i : pointIndex )
        totalArea += grid.areas.at(i);

      // iterate over values
      for(UInt i : pointIndex)
      {
        Double factor = multiplyWithArea ? std::pow(grid.points.at(i).r(), 2) : 1.0/totalArea;

        // compute area mean
        Vector data(grid.values.size());
        for(UInt j = 0; j < data.size(); j++)
          data(j) = grid.values.at(j).at(i);

        axpy(grid.areas.at(i) * factor, data, values.back());

        // compute RMS
        if(computeRms)
        {
          Vector tmpRms(grid.values.size());
          for(UInt j = 0; j < tmpRms.size(); j++)
            tmpRms(j) = std::pow(grid.values.at(j).at(i), 2);

          axpy(grid.areas.at(i) * factor, tmpRms, rms.back());
        }
      }
    }
    logTimerLoopEnd(times.size());

    // write file
    // ----------
    const UInt countDataFields = computeRms ? 2*maxDataFields+1 : maxDataFields + 1;
    Matrix A(values.size(), countDataFields);
    for(UInt k = 0; k<values.size(); k++)
    {
      A(k, 0) = gridTimes.at(k).mjd();
      copy(values.at(k).trans(), A.slice(k, 1, 1, values.at(k).size()));
    }
    if(removeMean)
    {
      logStatus<<"remove mean of time series"<<Log::endl;
      Vector meanValue(maxDataFields);
      for(UInt k = 0; k<values.size(); k++)
        axpy(1.0/static_cast<Double>(values.size()), A.slice(k, 1, 1, maxDataFields).trans(), meanValue);

      for(UInt k = 0; k<values.size(); k++)
        axpy(-1.0,  meanValue.trans(), A.slice(k, 1, 1, maxDataFields));
    }
    if(computeRms)
    {
      for(UInt k = 0; k < rms.size(); k++)
        for(UInt j = 0; j<rms.at(k).size(); j++)
          A(k, maxDataFields+1+j) = std::sqrt(rms.at(k)(j));
    }

    logStatus<<"write time series to file <"<<outputName<<">"<<Log::endl;
    InstrumentFile::write(outputName, Arc(A));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
