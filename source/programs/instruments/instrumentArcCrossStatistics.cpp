/***********************************************/
/**
* @file instrumentArcCrossStatistics.cpp
*
* @brief Compute RMS of an instrument time series or differences.
*
* @author Torsten Mayer-Guerr
* @date 2019-01-30
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Computes statistics of selected data columns between two \file{instrument files}{instrument} arc wise.
The \configFile{outputfileStatisticsTimeSeries}{instrument} contains for every arc one (mid) epoch
with statistics column(s). Possible statistics are
\begin{itemize}
  \item Correlation
  \begin{equation}
    \rho = \frac{\sum_i x_i y_i}{\sqrt{(\sum_i x_i^2) (\sum_i y_i^2})},
  \end{equation}
  \item Error RMS
  \begin{equation}
    rms = \sqrt{\frac{1}{N}\sum_i (x_i-y_i)^2},
  \end{equation}
  \item Nash-Sutcliffe coefficient (NSC)
  \begin{equation}
    nsc = 1- \frac{\sum_i (x_i-y_i)^2}{\sum_i (y_i-\bar{y})^2}.
  \end{equation}
\end{itemize}
With \config{removeArcMean} the mean of each data column of each arc is reduced before.

With \config{perColumn} separate statistics for each selected data column are computed,
otherwise an overall value is computed.

See also \program{InstrumentArcStatistics}, \program{InstrumentStatisticsTimeSeries}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Compute statistics for arcs of instrument data.
* @ingroup programsGroup */
class InstrumentArcCrossStatistics
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(InstrumentArcCrossStatistics, PARALLEL, "Compute RMS of an instrument time series or differences", Instrument, Statistics)

/***********************************************/

void InstrumentArcCrossStatistics::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName    fileNameInstrumentOut;
    FileName    fileNameInstrumentIn1, fileNameInstrumentIn2;
    UInt        startData, countData = MAX_UINT;
    Bool        removeMean, perColumn;
    std::string choice;
    std::function<double(const_MatrixSliceRef A, const_MatrixSliceRef B)> func;

    readConfig(config, "outputfileStatisticsTimeSeries", fileNameInstrumentOut, Config::MUSTSET, "", "statistics column(s) per arc, MISCVALUES");
    readConfig(config, "inputfileInstrument",            fileNameInstrumentIn1, Config::MUSTSET, "",   "");
    readConfig(config, "inputfileInstrumentReference",   fileNameInstrumentIn2, Config::MUSTSET, "",   "");
    if(readConfigChoice(config, "statistics", choice, Config::MUSTSET, "", ""))
    {
      if(readConfigChoiceElement(config, "correlation",   choice, ""))                                func = [] (const_MatrixSliceRef A, const_MatrixSliceRef B) {return inner(A, B)/rootMeanSquare(A)/rootMeanSquare(B);};
      if(readConfigChoiceElement(config, "errorRMS",      choice, "rms of differences"))              func = [] (const_MatrixSliceRef A, const_MatrixSliceRef B) {return rootMeanSquare(A-B);};
      if(readConfigChoiceElement(config, "nashSutcliffe", choice, "with respect to reference field")) func = [] (const_MatrixSliceRef A, const_MatrixSliceRef B) {return 1-quadsum(A-B)/quadsum(A-mean(A));};
      endChoice(config);
    }
    readConfig(config, "removeArcMean",   removeMean, Config::DEFAULT,  "1", "");
    readConfig(config, "startDataFields", startData,  Config::DEFAULT,  "0", "start");
    readConfig(config, "countDataFields", countData,  Config::OPTIONAL, "",  "number of data fields (default: all)");
    readConfig(config, "perColumn",       perColumn,  Config::DEFAULT,  "1", "compute statistic per column");
    if(isCreateSchema(config)) return;

    // =================================================

    logStatus<<"read instrument file <"<<fileNameInstrumentIn1<<">"<<Log::endl;
    logStatus<<"read instrument file <"<<fileNameInstrumentIn2<<">"<<Log::endl;
    InstrumentFile instrumentFile1(fileNameInstrumentIn1);
    InstrumentFile instrumentFile2(fileNameInstrumentIn2);
    InstrumentFile::checkArcCount({instrumentFile1, instrumentFile2});
    std::vector<Matrix> rows(instrumentFile1.arcCount());
    Parallel::forEach(rows, [&](UInt arcNo)
    {
      Arc arc1 = instrumentFile1.readArc(arcNo);
      Arc arc2 = instrumentFile2.readArc(arcNo);
      Arc::checkSynchronized({arc1, arc2});
      if(arc1.size() == 0)
        return Matrix();
      const Matrix A = arc1.matrix();
      const Matrix B = arc2.matrix();

      // compute statistics
      const UInt dataColumnCount = std::min(A.columns()-1-startData, countData);
      const UInt columnCount     = perColumn ? dataColumnCount : 1;
      Matrix row(1, 1+columnCount);
      row(0, 0) = (0.5*(arc1.front().time+arc1.back().time)+medianSampling(arc1.times())).mjd();
      for(UInt i=0; i<columnCount; i++)
      {
        const_MatrixSliceRef sliceA(A.column(1+startData+i, perColumn ? 1 : dataColumnCount));
        const_MatrixSliceRef sliceB(B.column(1+startData+i, perColumn ? 1 : dataColumnCount));
        row(0, 1+i) = (removeMean) ? func(sliceA-mean(sliceA), sliceB-mean(sliceB)) : func(sliceA, sliceB);
      }
      return row;
    }, comm);

    if(Parallel::isMaster(comm))
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

