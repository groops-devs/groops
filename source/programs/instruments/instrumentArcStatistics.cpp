/***********************************************/
/**
* @file instrumentArcStatistics.cpp
*
* @brief Compute statistics for arcs of instrument data.
*
* @author Sebastian Strasser
* @author Matthias Ellmer
* @author Torsten Mayer-Guerr
* @date 2016-07-18
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Computes statistics of selected data columns of \configFile{inputfileInstrument}{instrument} arc wise.
The \configFile{outputfileStatisticsTimeSeries}{instrument} contains for every arc one (mid) epoch
with statistics column(s). Possible statistics are root mean square, standard deviation,
mean, median, min, and max.

With \config{perColumn} separate statistics for each selected data column are computed,
otherwise an overall value is computed.

See also \program{InstrumentArcCrossStatistics}, \program{InstrumentStatisticsTimeSeries}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Compute statistics for arcs of instrument data.
* @ingroup programsGroup */
class InstrumentArcStatistics
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(InstrumentArcStatistics, PARALLEL, "Compute statistics for arcs of instrument data", Instrument, TimeSeries, Statistics)

/***********************************************/

void InstrumentArcStatistics::run(Config &config)
{
  try
  {
    FileName    fileNameInstrumentOut, fileNameInstrumentIn;
    UInt        startData, countData = MAX_UINT;
    Bool        perColumn, ignoreNan;
    std::string choice;
    std::function<double(const_MatrixSliceRef A)> func;

    readConfig(config, "outputfileStatisticsTimeSeries", fileNameInstrumentOut, Config::MUSTSET, "", "columns: mjd, statistics column(s) per instrument file");
    readConfig(config, "inputfileInstrument",            fileNameInstrumentIn,  Config::MUSTSET, "", "");
    if(readConfigChoice(config, "statistics", choice, Config::MUSTSET, "", ""))
    {
      if(readConfigChoiceElement(config, "rootMeanSquare",    choice)) func = rootMeanSquare;
      if(readConfigChoiceElement(config, "standardDeviation", choice)) func = standardDeviation;
      if(readConfigChoiceElement(config, "mean",              choice)) func = mean;
      if(readConfigChoiceElement(config, "median",            choice)) func = median;
      if(readConfigChoiceElement(config, "min",               choice)) func = static_cast<Double(*)(const_MatrixSliceRef)>(&min); // min and max are overloads, and the standard forbids automatic overload resolution in this context. So a manual cast is made in which the overload can be resolved successfully
      if(readConfigChoiceElement(config, "max",               choice)) func = static_cast<Double(*)(const_MatrixSliceRef)>(&max);
      if(readConfigChoiceElement(config, "epochCount",        choice)) func = [] (const_MatrixSliceRef A) {return static_cast<Double>(A.rows());};
      endChoice(config);
    }
    readConfig(config, "startDataFields", startData, Config::DEFAULT,  "0", "start");
    readConfig(config, "countDataFields", countData, Config::OPTIONAL, "",  "number of data fields (default: all)");
    readConfig(config, "perColumn",       perColumn, Config::DEFAULT,  "1", "compute statistic per column");
    readConfig(config, "ignoreNan",       ignoreNan, Config::DEFAULT,  "0", "ignore NaN values in input");
    if(isCreateSchema(config)) return;

    // =================================================

    logStatus<<"read instrument file <"<<fileNameInstrumentIn<<">"<<Log::endl;
    InstrumentFile instrumentFile(fileNameInstrumentIn);
    std::vector<Matrix> rows(instrumentFile.arcCount());
    Parallel::forEach(rows, [&](UInt arcNo)
    {
      const Arc    arc = instrumentFile.readArc(arcNo);
      const Matrix A   = arc.matrix();
      if(!A.size())
        return Matrix();

      // lambda function
      auto removeNan = [](const_MatrixSliceRef A)
      {
        std::vector<Double> data = flatten(A);
        data.erase(std::remove_if(data.begin(), data.end(), [](Double d) {return std::isnan(d);}), data.end());
        return Vector(data);
      };

      // compute statistics
      const UInt dataColumnCount = std::min(A.columns()-1-startData, countData);
      const UInt columnCount     = perColumn ? dataColumnCount : 1;
      Matrix row(1, 1+columnCount);
      row(0, 0) =  (0.5*(arc.front().time+arc.back().time+medianSampling(arc.times()))).mjd();
      for(UInt i=0; i<columnCount; i++)
      {
        const_MatrixSliceRef slice(A.column(1+startData+i, perColumn ? 1 : dataColumnCount));
        row(0, 1+i) = func(ignoreNan ? removeNan(slice) : slice);
      }
      return row;
    });

    if(Parallel::isMaster())
    {
      // copy rows to matrix
      rows.erase(std::remove_if(rows.begin(), rows.end(), [](const Matrix &A) {return !A.size();}), rows.end());
      Matrix A(rows.size(), rows.at(0).columns());
      for(UInt i=0; i<A.rows(); i++)
        copy(rows.at(i), A.row(i));

      logStatus<<"write arc statistics to instrument file <"<<fileNameInstrumentOut<<">"<<Log::endl;
      InstrumentFile::write(fileNameInstrumentOut, Arc(A));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
