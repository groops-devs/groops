/***********************************************/
/**
* @file instrumentStatisticsTimeSeries.cpp
*
* @brief Compute time series of statistics for intervals of instrument data.
*
* @author Sebastian Strasser
* @author Matthias Ellmer
* @date 2016-07-18
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring =
R"(
This program computes a time series of statistics for one or more instrument files.
Possible statistics are root mean square, standard deviation, mean, median, min, and max.
The columns of the output time series are defined either as one per \configFile{inputfileInstrument}{instrument}
or, if \config{perColumn} is true, statistics are computed per column for each file.
Providing e.g. 32 orbit files of GPS satellites results in a time series matrix
with columns: mjd, statisticsG01, statisticsG02, ..., statisticsG32.
If \config{intervals} are provided, the input data is split into these intervals
and one statistic is computed per interval. Otherwise, overall statistics are computed.
The instrument data considered for computation of the component-wise statistics
can be set with \config{startDataFields} and \config{countDataFields}.
The \config{factor} can be set to e.g. sqrt(3) to get 3D instead of 1D RMS values.

See also \program{InstrumentArcStatistics}, \program{InstrumentArcCrossStatistics}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"
#include "classes/timeSeries/timeSeries.h"

/***** CLASS ***********************************/

/** @brief Compute time series of statistics for intervals of instrument data.
* @ingroup programsGroup */
class InstrumentStatisticsTimeSeries
{
  static Double nanFunc(std::function<double(const_MatrixSliceRef A)> func, const_MatrixSliceRef A);

public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(InstrumentStatisticsTimeSeries, SINGLEPROCESS, "Compute time series of statistics for intervals of instrument data.", Instrument, TimeSeries, Statistics)

/***********************************************/

Double InstrumentStatisticsTimeSeries::nanFunc(std::function<double(const_MatrixSliceRef A)> func, const_MatrixSliceRef A)
{
  try
  {
    std::vector<Double> data = flatten(A);
    data.erase(std::remove_if(data.begin(), data.end(), [](Double d){ return std::isnan(d); }), data.end());
    return func(Vector(data));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void InstrumentStatisticsTimeSeries::run(Config &config)
{
  try
  {
    FileName              fileNameOutStatisticsTimeSeries;
    std::vector<FileName> fileNameInInstrument;
    TimeSeriesPtr         intervalsPtr;
    UInt                  startData = 0, countData = MAX_UINT;
    Bool                  perColumn, ignoreNan;
    Double                factor;
    std::string           choice;
    std::function<double(const_MatrixSliceRef A)> func;

    readConfig(config, "outputfileStatisticsTimeSeries", fileNameOutStatisticsTimeSeries, Config::MUSTSET, "", "columns: mjd, statistics column(s) per instrument file");
    readConfig(config, "inputfileInstrument",            fileNameInInstrument,            Config::MUSTSET, "", "");
    if(readConfigChoice(config, "statistics", choice, Config::MUSTSET, "", ""))
    {
      if(readConfigChoiceElement(config, "rootMeanSquare",    choice)) func = rootMeanSquare;
      if(readConfigChoiceElement(config, "standardDeviation", choice)) func = standardDeviation;
      if(readConfigChoiceElement(config, "mean",              choice)) func = mean;
      if(readConfigChoiceElement(config, "median",            choice)) func = median;
      if(readConfigChoiceElement(config, "sum",               choice)) func = sum;
      if(readConfigChoiceElement(config, "min",               choice)) func = static_cast<Double(*)(const_MatrixSliceRef)>(&min); // min and max are overloads, and the standard forbids automatic overload resolution in this context. So a manual cast is made in which the overload can be resolved successfully
      if(readConfigChoiceElement(config, "max",               choice)) func = static_cast<Double(*)(const_MatrixSliceRef)>(&max);
      if(readConfigChoiceElement(config, "epochCount",        choice)) func = [] (const_MatrixSliceRef A) { return static_cast<Double>(A.rows());};
      endChoice(config);
    }
    readConfig(config, "startDataFields", startData,    Config::DEFAULT,  "0", "start");
    readConfig(config, "countDataFields", countData,    Config::OPTIONAL, "",  "number of data fields (default: all)");
    readConfig(config, "perColumn",       perColumn,    Config::DEFAULT,  "0", "compute statistic per column");
    readConfig(config, "ignoreNan",       ignoreNan,    Config::DEFAULT,  "0", "ignore NaN values in statistic computation");
    readConfig(config, "intervals",       intervalsPtr, Config::DEFAULT,  "",  "intervals for statistics computation (one statistic per interval)");
    readConfig(config, "factor",          factor,       Config::DEFAULT,  "1", "e.g. sqrt(3) for 3D RMS");
    if(isCreateSchema(config)) return;

    std::vector<Vector> columns;

    // init time intervals
    std::vector<Time> intervals = intervalsPtr->times();
    if(intervals.size() >= 2)
    {
      columns.push_back(Vector(intervals.size()-1));
      for(UInt idInterval = 0; idInterval < intervals.size()-1; idInterval++)
        columns.back()(idInterval) = (0.5*(intervals.at(idInterval)+intervals.at(idInterval+1))).mjd();
    }
    else
      columns.push_back(Vector(1));

    // read instrument files and compute statistics
    for(UInt idInstrument = 0; idInstrument < fileNameInInstrument.size(); idInstrument++)
    {
      Arc arc;
      try
      {
        logStatus<<"read instrument file <"<<fileNameInInstrument.at(idInstrument)<<">"<<Log::endl;
        arc = InstrumentFile::read(fileNameInInstrument.at(idInstrument));
      }
      catch(std::exception &e)
      {
        logWarning<<e.what()<<" continue..."<<Log::endl;
        continue;
      }

      if(arc.size() == 0)
        continue;

      // determine data spans
      Matrix A = arc.matrix();
      std::vector<UInt> intervalEpochIds;
      if(intervals.size() >= 2)
      {
        std::vector<Double> times = Vector(A.column(0));
        for(UInt i = 0; i < intervals.size(); i++)
        {
          const auto it = std::find_if(times.begin(), times.end(), [intervals, i] (Double t) { return t >= intervals.at(i).mjd(); });
          intervalEpochIds.push_back(std::distance(times.begin(), it));
        }
      }
      else
        intervalEpochIds = {0, arc.size()};

      // compute statistics
      const UInt dataColumnCount = std::min(A.columns()-1-startData, countData);
      const UInt rowCount        = intervalEpochIds.size()-1;
      const UInt columnCount     = perColumn ? dataColumnCount : 1;
      for(UInt idCol = 0; idCol < columnCount; idCol++)
      {
        Vector column(rowCount);
        for(UInt idRow = 0; idRow < rowCount; idRow++)
        {
          const UInt rows = intervalEpochIds.at(idRow+1)-intervalEpochIds.at(idRow);
          if(rows)
          {
            Matrix slice = A.slice(intervalEpochIds.at(idRow), 1+startData+idCol, rows, perColumn ? 1 : dataColumnCount);
            column(idRow) = factor * (ignoreNan ? nanFunc(func, slice) : func(slice));
          }
          else
            column(idRow) = NAN_EXPR;
        }
        columns.push_back(column);
      }
    }

    // save file
    logStatus<<"write statistics time series to file <"<<fileNameOutStatisticsTimeSeries<<">"<<Log::endl;
    Matrix A(columns.at(0).size(), columns.size());
    for(UInt i = 0; i < columns.size(); i++)
      copy(columns.at(i), A.column(i));
    writeFileMatrix(fileNameOutStatisticsTimeSeries, A);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
