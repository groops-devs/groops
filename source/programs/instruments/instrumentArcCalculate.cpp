/***********************************************/
/**
* @file instrumentArcCalculate.cpp
*
* @brief Instrument data calculation arc wise.
*
* @author Torsten Mayer-Guerr
* @date 2018-05-01
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program manipulates the data columns every arc of an \file{instrument file}{instrument} similar to
\program{FunctionsCalculate}, see there for more details.
If several \configFile{inputfileInstrument}{instrument}s are given the data columns are copied side by side.
For this the instrument files must be synchronized (see \program{InstrumentSynchronize}). For the data
columns the standard data variables are available, see~\reference{dataVariables}{general.parser:dataVariables}.
For the time column (MJD) a variable \verb|epoch| (together with \verb|epochmean|, \verb|epochmin|, \ldots)
is defined additionally.

The content of \configFile{outputfileInstrument}{instrument} is controlled by \config{outColumn}.
The number of \config{outColumn} must agree with the selected \configClass{outType}{instrumentTypeType}.
The algorithm to compute the output is as follows:
The expressions in \config{outColumn} are evaluated once for each epoch of the input.
The variables \verb|data0|,~\verb|data1|,~\ldots are replaced by the according values from the input columns before.
If no \config{outColumn} are specified all input columns are used instead directly.
The \configClass{instrument type}{instrumentTypeType} can be specified with \config{outType} and must be agree with the number of columns.

An extra \config{statistics} file can be generated with one mid epoch per arc. For the computation of the \config{outColumn} values
all~\reference{dataVariables}{general.parser:dataVariables} are available (e.g. \verb|epochmin|, \verb|data0mean|, \verb|data1std|, \ldots)
inclusively the \config{constant}s and estimated \config{parameter}s but without the \verb|data0|,~\verb|data1|,~\ldots itself.
The variables and the numbering of the columns refers to the \configFile{outputfileInstrument}{instrument}.

See also \program{FunctionsCalculate}, \program{MatrixCalculate}.
)";

/***********************************************/

#include "programs/program.h"
#include "parser/dataVariables.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"

/***********************************************/

/** @brief Instrument data calculation arc wise.
* @ingroup programsGroup */
class InstrumentArcCalculate
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(InstrumentArcCalculate, PARALLEL, "Instrument data calculation arc wise.", Instrument)

/***********************************************/

void InstrumentArcCalculate::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName                           fileNameOut, fileNameStatistics;
    std::vector<FileName>              fileNamesIn;
    std::vector<ExpressionVariablePtr> constExpr, paramExpr;
    std::vector<ExpressionVariablePtr> lsaExpr, removeExpr;
    std::vector<ExpressionVariablePtr> outExpr, statisticsExpr;
    Epoch::Type                        type = Epoch::EMPTY;

    readConfig(config, "outputfileInstrument",    fileNameOut, Config::OPTIONAL, "",      "");
    readConfig(config, "inputfileInstrument",     fileNamesIn, Config::MUSTSET,  "",      "data columns are appended to the right");
    readConfig(config, "constant",                constExpr,   Config::OPTIONAL, "",      "define a constant by name=value");
    readConfig(config, "parameter",               paramExpr,   Config::OPTIONAL, "",      "define a parameter by name[=value]");
    readConfig(config, "leastSquares",            lsaExpr,     Config::OPTIONAL, "",      "try to minimize the expression by adjustment of the parameters");
    readConfig(config, "removalCriteria",         removeExpr,  Config::OPTIONAL, "",      "row is removed if one criterion evaluates true.");
    readConfig(config, "outType",                 type,        Config::OPTIONAL, "",      "");
    readConfig(config, "outColumn",               outExpr,     Config::OPTIONAL, "data0", "expression of output columns, extra 'epoch' variable");
    if(readConfigSequence(config, "statistics", Config::OPTIONAL, "", ""))
    {
      readConfig(config, "outputfileInstrument", fileNameStatistics, Config::MUSTSET, "",         "instrument file with mid epoch per arc, data columns are user defined");
      readConfig(config, "outColumn",            statisticsExpr,     Config::MUSTSET, "data0rms", "expression to compute statistics columns, data* are from outColumn");
      endSequence(config);
    }
    if(isCreateSchema(config)) return;

    // open files
    // ----------
    std::vector<InstrumentFile> file(fileNamesIn.size());
    for(UInt i=0; i<file.size(); i++)
      file.at(i).open(fileNamesIn.at(i));
    for(UInt i=1; i<file.size(); i++)
      InstrumentFile::checkArcCount({file.at(0), file.at(i)});

    // create data variables
    // ---------------------
    VariableList varListGlobal;
    // get real variable names, otherwise all named after config element
    std::for_each(constExpr.begin(), constExpr.end(), [&](auto expr) {expr->parseVariableName();});
    std::for_each(paramExpr.begin(), paramExpr.end(), [&](auto expr) {expr->parseVariableName();});
    std::for_each(constExpr.begin(), constExpr.end(), [&](auto expr) {varListGlobal.addVariable(expr);});
    std::for_each(paramExpr.begin(), paramExpr.end(), [&](auto expr) {varListGlobal.addVariable(expr);});

    // arc wise computation
    // --------------------
    logStatus<<"computing arcs"<<Log::endl;
    std::vector<Arc> arcList(file.at(0).arcCount());
    Matrix statistics(file.at(0).arcCount(), 1+statisticsExpr.size());
    Parallel::forEach(arcList, [&](UInt arcNo)
    {
      // read data
      std::vector<Arc> arc(file.size());
      for(UInt i=0; i<arc.size(); i++)
        arc.at(i) = file.at(i).readArc(arcNo);
      for(UInt i=1; i<arc.size(); i++)
        Arc::checkSynchronized({arc.at(0), arc.at(i)});

      // copy data to one matrix + extra time vector
      std::vector<Time> times = (arc.at(0).times());
      std::vector<UInt> index(1, 0);
      for(UInt i=0; i<arc.size(); i++)
        index.push_back( arc.at(i).at(0).data().rows() + index.back() );
      Matrix data(arc.at(0).size(), index.back());
      for(UInt i=0; i<arc.size(); i++)
        copy(arc.at(i).matrix().column(1, index.at(i+1)-index.at(i)), data.column(index.at(i), index.at(i+1)-index.at(i))); // without time column

      auto varListArc       = varListGlobal;
      auto varListArcWoData = varListGlobal;
      addDataVariables("epoch", times, varListArc);
      addDataVariables(data, varListArc);

      // least squares adjustment
      if(lsaExpr.size())
      {
        Vector l(data.rows()*lsaExpr.size());
        Matrix A(data.rows()*lsaExpr.size(), paramExpr.size());
        for(UInt k=0; k<lsaExpr.size(); k++)
        {
          std::vector<ExpressionVariablePtr> designExpr(paramExpr.size());
          for(UInt s=0; s<paramExpr.size(); s++)
          {
            designExpr.at(s) = lsaExpr.at(k)->derivative(paramExpr.at(s)->name(), varListArc);
            designExpr.at(s)->simplify(varListArc);
          }

          for(UInt i=0; i<data.rows(); i++)
          {
            evaluateDataVariables(data, i, varListArc);
            varListArc.setVariable("epoch", times.at(i).mjd());
            l(i+k*data.rows()) = -lsaExpr.at(k)->evaluate(varListArc); // observations
            for(UInt s=0; s<designExpr.size(); s++)
              A(i+k*data.rows(),s) = designExpr.at(s)->evaluate(varListArc); // columns of design matrix
          }
          undefineDataVariables(data, varListArc);
          varListArc.undefineVariable("epoch");
        }

        Vector x = leastSquares(A, l);
        for(UInt s=0; s<paramExpr.size(); s++)
        {
          x(s) += paramExpr.at(s)->evaluate(varListArc);
          varListArc.setVariable(paramExpr.at(s)->name(),  x(s) );
          varListArcWoData.setVariable(paramExpr.at(s)->name(),  x(s) );
        }
      } // if(lsa)

      // create output arc
      std::vector<Time> timesOut(data.rows());
      Matrix outData(data.rows(), 1 + (outExpr.size() ? outExpr.size() : data.columns())); // first column for time
      UInt row = 0;
      for(UInt i=0; i<outData.rows(); i++)
      {
        evaluateDataVariables(data, i, varListArc);
        varListArc.setVariable("epoch", times.at(i).mjd());
        if(!std::any_of(removeExpr.begin(), removeExpr.end(), [&](auto expr) {return expr->evaluate(varListArc) != 0;}))
        {
          timesOut.at(row) = times.at(i);
          for(UInt k=0; k<outData.columns()-1; k++)
            outData(row, 1+k) = outExpr.size() ? outExpr.at(k)->evaluate(varListArc) : data(i, k);
          row++;
        }
      }
      timesOut.resize(row);
      if(row < outData.rows())
        outData = outData.row(0, row);

      // arc statistics
      if(statisticsExpr.size())
      {
        auto varList = varListArcWoData;
        addDataVariables("epoch", timesOut, varList);
        addDataVariables(outData.column(1, outData.columns()-1), varList);
        statistics(arcNo, 0) = (0.5*(timesOut.front()+timesOut.back()+medianSampling(timesOut))).mjd();
        for(UInt k=0; k<statisticsExpr.size(); k++)
          statistics(arcNo, 1+k) = statisticsExpr.at(k)->evaluate(varList);
      }

      return Arc(timesOut, outData, type);
    }, comm);

    // write results
    // -------------
    if(!fileNameOut.empty() && Parallel::isMaster(comm))
    {
      logStatus<<"write instrument data <"<<fileNameOut<<">"<<Log::endl;
      InstrumentFile::write(fileNameOut, arcList);
      Arc::printStatistics(arcList);
    }

    if(!fileNameStatistics.empty() && statistics.size())
    {
      Parallel::reduceSum(statistics, 0, comm);
      if(Parallel::isMaster(comm))
      {
        logStatus<<"write statistics file <"<<fileNameStatistics<<">"<<Log::endl;
        InstrumentFile::write(fileNameStatistics, Arc(statistics));
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
