/***********************************************/
/**
* @file griddedDataCalculate.cpp
*
* @brief Calculate values of gridded data.
*
* @author Torsten Mayer-Guerr
* @date 2012-11-20
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program manipulates \file{grid files}{griddedData} with data in columns similar to
\program{FunctionsCalculate}, see there for more details.
If several \config{inputfile}s are given the data columns are copied side by side.
All \config{inputfile}s must contain the same grid points.
The columns are enumerated by \verb|data0|,~\verb|data1|,~\ldots.

The content of \configFile{outputfileGriddedData}{griddedData} is controlled by \config{outColumn}.
The algorithm to compute the output is as follows:
The expressions in \config{outColumn} are evaluated once for each grid point of the input.
The variables \verb|data0|,~\verb|data1|,~\ldots are replaced by the according values
from the input columns before.
Additional variables are available, e.g. \verb|index|, \verb|data0rms|,
see~\reference{dataVariables}{general.parser:dataVariables}.

For a simplified handling \config{constant}s can be defined by \verb|name=value|.
It is also possible to estimate \config{parameter}s in a least squares adjustment.
The \config{leastSquares} serves as template for observation equations for every point.
The expression \config{leastSquares} is evaluated for each grid point.
The variables \verb|data0|,~\verb|data1|,~\ldots are replaced by the according values from the input columns before.
In the next step the parameters are estimated in order to minimize the expressions in \config{leastSquares}
in the sense of least squares.

Afterwards grid points are removed if one of the \config{removalCriteria} expressions
for this grid point evaluates true (not zero).

An extra \configFile{statistics:outputfile}{matrix} can be generated with one row of data.
For the computation of the \config{outColumn} values
all~\reference{dataVariables}{general.parser:dataVariables} are available
(e.g. \verb|data3mean|, \verb|data4std|) inclusively the \config{constant}s and
estimated \config{parameter}s but without the \verb|data0|,~\verb|data1|,~\ldots itself.
The variables and the numbering of the columns refers to the \configFile{outputfileGriddedData}{griddedData}.

See also \program{FunctionsCalculate}, \program{InstrumentArcCalculate}, \program{MatrixCalculate}.
)";

/***********************************************/

#include "programs/program.h"
#include "parser/dataVariables.h"
#include "files/fileMatrix.h"
#include "files/fileGriddedData.h"
#include "misc/miscGriddedData.h"

/***********************************************/

/** @brief Calculate values of gridded data.
* @ingroup programsGroup */
class GriddedDataCalculate
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GriddedDataCalculate, SINGLEPROCESS, "Calculate values of gridded data.", Grid)

/***********************************************/

void GriddedDataCalculate::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName                           fileNameOut, fileNameStatistics;
    std::vector<FileName>              fileNamesIn;
    std::vector<ExpressionVariablePtr> constExpr, paramExpr;
    std::vector<ExpressionVariablePtr> lsaExpr, removeExpr;
    std::vector<ExpressionVariablePtr> valueExpr, statisticsExpr;
    ExpressionVariablePtr              lonExpr, latExpr, heightExpr, areaExpr;
    Bool                               computeArea;
    Double                             a, f;

    readConfig(config, "outputfileGriddedData", fileNameOut, Config::OPTIONAL, "", "");
    readConfig(config, "inputfileGriddedData",  fileNamesIn, Config::MUSTSET,  "", "");
    readConfig(config, "constant",              constExpr,   Config::OPTIONAL, "", "define a constant by name=value");
    readConfig(config, "parameter",             paramExpr,   Config::OPTIONAL, "", "define a parameter by name[=value]");
    readConfig(config, "leastSquares",          lsaExpr,     Config::OPTIONAL, "", "try to minimize the expression by adjustment of the parameters");
    readConfig(config, "removalCriteria",       removeExpr,  Config::OPTIONAL, "", "points are removed if one criterion evaluates true. data0 is the first data field.");
    readConfig(config, "longitude",             lonExpr,     Config::MUSTSET,  "longitude", "expression");
    readConfig(config, "latitude",              latExpr,     Config::MUSTSET,  "latitude",  "expression");
    readConfig(config, "height",                heightExpr,  Config::MUSTSET,  "height",    "expression");
    readConfig(config, "area",                  areaExpr,    Config::OPTIONAL, "area",      "expression: e.g. deltaL * 2.0 * sin(deltaB/2.0) * cos(latitude/rho)");
    readConfig(config, "value",                 valueExpr,   Config::OPTIONAL, "data0",     "expression to compute values (input columns are named data0, data1, ...)");
    readConfig(config, "computeArea",           computeArea, Config::DEFAULT,  "0",         "automatically area computation of rectangular grids (overwrite area)");
    readConfig(config, "R",                     a,           Config::DEFAULT,  STRING_DEFAULT_GRS80_a, "reference radius for ellipsoidal coordinates");
    readConfig(config, "inverseFlattening",     f,           Config::DEFAULT,  STRING_DEFAULT_GRS80_f, "reference flattening for ellipsoidal coordinates");
    if(readConfigSequence(config, "statistics", Config::OPTIONAL, "", ""))
    {
      readConfig(config, "outputfile", fileNameStatistics, Config::MUSTSET, "",         "matrix with one row, columns are user defined");
      readConfig(config, "outColumn",  statisticsExpr,     Config::MUSTSET, "data0rms", "expression to compute statistics columns, data* are the outputColumns");
      endSequence(config);
    }
    if(isCreateSchema(config)) return;

    // =====================================================

    // read data
    // ---------
    GriddedData gridIn;
    for(UInt idFile=0; idFile<fileNamesIn.size(); idFile++)
    {
      logStatus<<"reading grid from file <"<<fileNamesIn.at(idFile)<<">"<<Log::endl;
      GriddedData grid;
      readFileGriddedData(fileNamesIn.at(idFile), grid);

      // test grid equivalence
      if(gridIn.points.size() != 0)
      {
        if(gridIn.points.size() != grid.points.size())
          throw(Exception("grid differ in point count"));
        for(UInt i=0; i<gridIn.points.size(); i++)
          if(std::fabs((gridIn.points.at(i)-grid.points.at(i)).r()/gridIn.points.at(i).r())>1e-10)
          {
            logWarning<<"grid points differ: ignored"<<Log::endl;
            break;
          }
        // append values
        gridIn.values.insert(gridIn.values.end(), grid.values.begin(), grid.values.end());
      }
      else
        gridIn = grid;
    } // for(idFile)

    // =====================================================

    // create data variables
    // ---------------------
    auto varList = config.getVarList();
    // get real variable names, otherwise all named after config element
    std::for_each(constExpr.begin(), constExpr.end(), [&](auto expr) {expr->parseVariableName(varList);});
    std::for_each(paramExpr.begin(), paramExpr.end(), [&](auto expr) {expr->parseVariableName(varList);});

    std::set<std::string> usedVariables;
    lonExpr   ->usedVariables(varList, usedVariables);
    latExpr   ->usedVariables(varList, usedVariables);
    heightExpr->usedVariables(varList, usedVariables);
    if(areaExpr)
      areaExpr->usedVariables(varList, usedVariables);
    std::for_each(valueExpr.begin(),  valueExpr.end(),  [&](auto expr) {expr->usedVariables(varList, usedVariables);});
    std::for_each(lsaExpr.begin(),    lsaExpr.end(),    [&](auto expr) {expr->usedVariables(varList, usedVariables);});
    std::for_each(removeExpr.begin(), removeExpr.end(), [&](auto expr) {expr->usedVariables(varList, usedVariables);});
    std::for_each(constExpr.begin(), constExpr.end(), [&](auto expr) {addVariable(expr, varList);});
    std::for_each(paramExpr.begin(), paramExpr.end(), [&](auto expr) {addVariable(expr, varList);});
    auto varListWoData = varList;
    addDataVariables(gridIn, varList, usedVariables);

    // =====================================================

    // build observation equations
    // ---------------------------
    if(lsaExpr.size())
    {
      logStatus<<"least squares adjustment"<<Log::endl;
      Vector l(gridIn.points.size()*lsaExpr.size());
      Matrix A(gridIn.points.size()*lsaExpr.size(), paramExpr.size());

      for(UInt k=0; k<lsaExpr.size(); k++)
      {
        std::vector<ExpressionVariablePtr> lsaA(paramExpr.size());
        for(UInt s=0; s<paramExpr.size(); s++)
        {
          lsaA.at(s) = lsaExpr.at(k)->derivative(paramExpr.at(s)->name(), varList);
          lsaA.at(s)->simplify(varList);
        }
        lsaExpr.at(k)->simplify(varList);

        for(UInt i=0; i<gridIn.points.size(); i++)
        {
          evaluateDataVariables(gridIn, i, varList);
          l(i+k*gridIn.points.size()) = -lsaExpr.at(k)->evaluate(varList); // observations
          for(UInt s=0; s<lsaA.size(); s++)
            A(i+k*gridIn.points.size(),s) = lsaA.at(s)->evaluate(varList); // columns of design matrix
        }
        undefineDataVariables(gridIn, varList);
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

    // =====================================================

    // calculate output grid
    // ---------------------
    logStatus<<"calculate output matrix"<<Log::endl;
    lonExpr->simplify(varList);
    latExpr->simplify(varList);
    heightExpr->simplify(varList);
    if(areaExpr)
      areaExpr->simplify(varList);
    std::for_each(removeExpr.begin(), removeExpr.end(), [&](auto expr) {expr->simplify(varList);});
    std::for_each(valueExpr.begin(),  valueExpr.end(),  [&](auto expr) {expr->simplify(varList);});

    GriddedData gridOut;
    gridOut.ellipsoid = Ellipsoid(a,f);
    gridOut.values.resize(valueExpr.size());
    for(UInt i=0; i<gridIn.points.size(); i++)
    {
      evaluateDataVariables(gridIn, i, varList);
      if(std::any_of(removeExpr.begin(), removeExpr.end(), [&](auto expr) {return expr->evaluate(varList) != 0.;}))
        continue;
      Double L = lonExpr->evaluate(varList);
      Double B = latExpr->evaluate(varList);
      Double h = heightExpr->evaluate(varList);
      gridOut.points.push_back( gridOut.ellipsoid(Angle(L*DEG2RAD), Angle(B*DEG2RAD), h) );
      if(areaExpr)
        gridOut.areas.push_back( areaExpr->evaluate(varList) );
      for(UInt k=0; k<valueExpr.size(); k++)
        gridOut.values.at(k).push_back( valueExpr.at(k)->evaluate(varList) );
    }

    // =====================================================

    if(computeArea)
    {
      logStatus<<"compute area"<<Log::endl;
      if(!gridOut.computeArea())
        logWarning<<"Compute areas only possible with regular rectangular grids"<<Log::endl;
    }

    // =====================================================

    // save grid
    // ---------
    if(!fileNameOut.empty())
    {
      logStatus<<"save grid <"<<fileNameOut<<">"<<Log::endl;
      writeFileGriddedData(fileNameOut, gridOut);
      MiscGriddedData::printStatistics(gridOut);
    }

    // statistics
    // ----------
    if(!fileNameStatistics.empty())
    {
      logStatus<<"write statistics to <"<<fileNameStatistics<<">"<<Log::endl;
      auto varList = varListWoData;
      std::set<std::string> usedVariables;
      std::for_each(statisticsExpr.begin(), statisticsExpr.end(), [&](auto expr) {expr->usedVariables(varList, usedVariables);});
      addDataVariables(gridOut, varList, usedVariables);
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
