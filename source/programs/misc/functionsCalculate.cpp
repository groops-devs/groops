/***********************************************/
/**
* @file functionsCalculate.cpp
*
* @brief Functions Calculate.
*
* @author Torsten Mayer-Guerr
* @date 2009-09-11
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program manipulates \file{matrix files}{matrix} with data in columns.
If several \config{inputfile}s are given the data columns are copied side by side.
All \config{inputfile}s must contain the same number of rows.
The columns are enumerated by \verb|data0|,~\verb|data1|,~\ldots.

The content of \configFile{outputfile}{matrix} is controlled by \config{outColumn}.
The algorithm to compute the output is as follows:
The expressions in \config{outColumn} are evaluated once for each row of the input.
The variables \verb|data0|,~\verb|data1|,~\ldots are replaced by the according values from the input columns before.
Additional variables are available, e.g. \verb|index|, \verb|data0rms|, see~\reference{dataVariables}{general.parser:dataVariables}.
If no \config{outColumn} are specified all input columns are used instead directly.

For a simplified handling \config{constant}s can be defined by \verb|name=value|, e.g. \verb|annual=365.25|.
It is also possible to estimate \config{parameter}s in a least squares adjustment.
The \config{leastSquares} serves as template for observation equations for every row.
The expression \config{leastSquares} is evaluated for each row in the \config{inputfile}.
The variables \verb|data0|,~\verb|data1|,~\ldots are replaced by the according values from the input columns before.
In the next step the parameters are estimated in order to minimize the expressions in \config{leastSquares}
in the sense of least squares.

Afterwards complete rows are removed if one of the \config{removalCriteria} expressions for this row evaluates true (not zero).

An extra \config{statistics} file can be generated with one row of data. For the computation of the \config{outColumn} values
all~\reference{dataVariables}{general.parser:dataVariables} are available (e.g. \verb|data3mean|, \verb|data4std|)
inclusively the \config{constant}s and estimated \config{parameter}s but without the \verb|data0|,~\verb|data1|,~\ldots itself.
The variables and the numbering of the columns refers to the \configFile{outputfile}{matrix}.

First example: To calculate the mean of two values at each row set \config{outColumn} to \verb|0.5*(data1+data0)|.

Second example: An input file contain a column with times and a column with values.
To remove a trend from the values define the \config{parameter}s \verb|trend| and \verb|bias|.
The observation equation in \config{leastSquares} is \verb|data1 - (trend*data0+bias)|.
For output you can define the following columns for example:
\begin{itemize}
\item \config{outColumn}=\verb|data0|: points in time.
\item \config{outColumn}=\verb|data1|: the values itself.
\item \config{outColumn}=\verb|trend*data0+bias|: the linear fit.
\item \config{outColumn}=\verb|data1-trend*data0-bias|: the residuals.
\end{itemize}
The extra statistics file could contain in this case:
\begin{itemize}
\item \config{outColumn}=\verb|data0max-data0min|: time span.
\item \config{outColumn}=\verb|bias|: estimated parameter.
\item \config{outColumn}=\verb|trend|: estimated parameter.
\item \config{outColumn}=\verb|data3rms|: root mean square of the residuals.
\end{itemize}

See also \program{InstrumentArcCalculate}, \program{GriddedDataCalculate}, \program{MatrixCalculate}.
)";

/***********************************************/

#include "programs/program.h"
#include "parser/dataVariables.h"
#include "files/fileMatrix.h"

/***********************************************/

/** @brief Functions Calculate.
* @ingroup programsGroup */
class FunctionsCalculate
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(FunctionsCalculate, SINGLEPROCESS, "Functions Calulate", Misc, Matrix, TimeSeries)

/***********************************************/

void FunctionsCalculate::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName                           fileNameOut, fileNameStatistics;
    std::vector<FileName>              fileNamesIn;
    std::vector<ExpressionVariablePtr> constExpr, paramExpr;
    std::vector<ExpressionVariablePtr> lsaExpr, removeExpr;
    std::vector<ExpressionVariablePtr> outExpr, statisticsExpr;

    readConfig(config, "outputfile",      fileNameOut, Config::MUSTSET,  "",      "");
    readConfig(config, "inputfile",       fileNamesIn, Config::MUSTSET,  "",      "");
    readConfig(config, "constant",        constExpr,   Config::OPTIONAL, "",      "define a constant by name=value");
    readConfig(config, "parameter",       paramExpr,   Config::OPTIONAL, "",      "define a parameter by name[=value]");
    readConfig(config, "leastSquares",    lsaExpr,     Config::OPTIONAL, "",      "try to minimize the expression by adjustment of the parameters");
    readConfig(config, "removalCriteria", removeExpr,  Config::OPTIONAL, "",      "row is removed if one criterion evaluates true.");
    readConfig(config, "outColumn",       outExpr,     Config::OPTIONAL, "data0", "expression to compute output columns (input columns are named data0, data1, ...)");
    if(readConfigSequence(config, "statistics", Config::OPTIONAL, "", ""))
    {
      readConfig(config, "outputfile", fileNameStatistics, Config::MUSTSET, "",         "matrix with one row, columns are user defined");
      readConfig(config, "outColumn",  statisticsExpr,     Config::MUSTSET, "data0rms", "expression to compute statistics columns, data* are the outputColumns");
      endSequence(config);
    }
    if(isCreateSchema(config)) return;

    // read data
    // ---------
    Matrix data;
    {
      std::vector<Matrix> data2(fileNamesIn.size());
      for(UInt i=0; i<fileNamesIn.size(); i++)
      {
        logStatus<<"read input from <"<<fileNamesIn.at(i)<<">"<<Log::endl;
        readFileMatrix(fileNamesIn.at(i), data2.at(i));
      }
      // test data
      for(UInt i=1; i<data2.size(); i++)
        if(data2.at(i).rows() != data2.at(0).rows())
          throw(Exception("all input data must have the same count of rows"));
      UInt cols = 0;
      for(UInt i=0; i<data2.size(); i++)
        cols += data2.at(i).columns();
      data = Matrix(data2.at(0).rows(), cols);
      UInt idx = 0;
      for(UInt i=0; i<data2.size(); i++)
      {
        copy(data2.at(i), data.column(idx,data2.at(i).columns()));
        idx += data2.at(i).columns();
      }
    }

    // create data variables
    // ---------------------
    auto varList = config.getVarList();
    // get real variable names, otherwise all named after config element
    std::for_each(constExpr.begin(), constExpr.end(), [&](auto expr) {expr->parseVariableName(varList);});
    std::for_each(paramExpr.begin(), paramExpr.end(), [&](auto expr) {expr->parseVariableName(varList);});

    std::set<std::string> usedVariables;
    std::for_each(outExpr.begin(),    outExpr.end(),    [&](auto expr) {expr->usedVariables(varList, usedVariables);});
    std::for_each(lsaExpr.begin(),    lsaExpr.end(),    [&](auto expr) {expr->usedVariables(varList, usedVariables);});
    std::for_each(removeExpr.begin(), removeExpr.end(), [&](auto expr) {expr->usedVariables(varList, usedVariables);});
    std::for_each(constExpr.begin(), constExpr.end(), [&](auto expr) {addVariable(expr, varList);});
    std::for_each(paramExpr.begin(), paramExpr.end(), [&](auto expr) {addVariable(expr, varList);});
    auto varListWoData = varList;
    addDataVariables(data, varList, usedVariables);

    // build observation equations
    // ---------------------------
    if(lsaExpr.size())
    {
      logStatus<<"least squares adjustment"<<Log::endl;
      Vector l(data.rows()*lsaExpr.size());
      Matrix A(data.rows()*lsaExpr.size(), paramExpr.size());

      for(UInt k=0; k<lsaExpr.size(); k++)
      {
        std::vector<ExpressionVariablePtr> lsaA(paramExpr.size());
        for(UInt s=0; s<paramExpr.size(); s++)
        {
          lsaA.at(s) = lsaExpr.at(k)->derivative(paramExpr.at(s)->name(), varList);
          lsaA.at(s)->simplify(varList);
        }
        lsaExpr.at(k)->simplify(varList);

        for(UInt i=0; i<data.rows(); i++)
        {
          evaluateDataVariables(data, i, varList);
          l(i+k*data.rows()) = -lsaExpr.at(k)->evaluate(varList);  // observations
          for(UInt s=0; s<lsaA.size(); s++)
            A(i+k*data.rows(), s) = lsaA.at(s)->evaluate(varList);  // columns of design matrix
        }
        undefineDataVariables(data, varList);
      }

      Vector x = leastSquares(A,l);
      for(UInt s=0; s<paramExpr.size(); s++)
      {
        x(s) += paramExpr.at(s)->evaluate(varList);
        paramExpr.at(s)->setValue( x(s) );
        varList[paramExpr.at(s)->name()]->setValue( x(s) );
        varListWoData[paramExpr.at(s)->name()]->setValue( x(s) );
        logInfo<<"  "<<paramExpr.at(s)->name()<<" = "<<x(s)<<Log::endl;
      }
    }

    // calculate output matrix
    // -----------------------
    logStatus<<"calculate output matrix"<<Log::endl;
    std::for_each(outExpr.begin(),    outExpr.end(),    [&](auto expr) {expr->simplify(varList);});
    std::for_each(removeExpr.begin(), removeExpr.end(), [&](auto expr) {expr->simplify(varList);});
    Matrix outData(data.rows(), outExpr.size() ? outExpr.size() : data.columns());
    UInt row = 0;
    for(UInt i=0; i<outData.rows(); i++)
    {
      evaluateDataVariables(data, i, varList);
      if(!std::any_of(removeExpr.begin(), removeExpr.end(), [&](auto expr) {return expr->evaluate(varList) != 0;}))
      {
        for(UInt k=0; k<outData.columns(); k++)
          outData(row, k) = outExpr.size() ? outExpr.at(k)->evaluate(varList) : data(i, k);
        row++;
      }
    }
    if(row < outData.rows())
    {
      logInfo<<"  "<<outData.rows()-row<<" rows removed"<<Log::endl;
      outData = outData.row(0, row);
    }

    // write output
    // ------------
    if(!fileNameOut.empty())
    {
      logStatus<<"write output to <"<<fileNameOut<<">"<<Log::endl;
      writeFileMatrix(fileNameOut, outData);
    }

    // statistics
    // ----------
    if(!fileNameStatistics.empty())
    {
      logStatus<<"write statistics to <"<<fileNameStatistics<<">"<<Log::endl;
      auto varList = varListWoData;
      std::set<std::string> usedVariables;
      std::for_each(statisticsExpr.begin(), statisticsExpr.end(), [&](auto expr) {expr->usedVariables(varList, usedVariables);});
      addDataVariables(outData, varList, usedVariables);
      Matrix statistics(1, statisticsExpr.size());
      for(UInt k=0; k<statistics.columns(); k++)
        statistics(0, k) = statisticsExpr.at(k)->evaluate(varList);
      writeFileMatrix(fileNameStatistics, statistics);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
