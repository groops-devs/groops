/***********************************************/
/**
* @file griddedDataTimeSeries2PotentialCoefficients.cpp
*
* @brief Convert a griddedDataTimeSeries to a sequence of PotentialCoefficients files.
*
* @author Torsten Mayer-Guerr
* @date 2023-11-08
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program estimate potential coefficients from
\configFile{inputfileGriddedDataTimeSeries}{griddedDataTimeSeries}
in the same way as \program{GriddedData2PotentialCoefficients}
but not only for one grid but for each epoch of
\configClass{timeSeries}{timeSeriesType} of if not set
for the temporal nodal points from the inputfile.
The \configFile{outputfilePotentialCoefficients}{potentialCoefficients}
(one for each \config{value}) are written for each epoch with the expansion
of \config{variableLoopTime} and \config{variableLoopIndex}
(see \reference{text parser}{general.parser:text}).

See also \program{GriddedData2PotentialCoefficients}.
)";

/***********************************************/

#include "programs/program.h"
#include "parser/dataVariables.h"
#include "files/fileGriddedDataTimeSeries.h"
#include "files/fileSphericalHarmonics.h"
#include "classes/timeSeries/timeSeries.h"
#include "classes/kernel/kernel.h"
#include "misc/miscGriddedData.h"

/***** CLASS ***********************************/

/** @brief Convert a griddedDataTimeSeries to a sequence of PotentialCoefficients files.
* @ingroup programsGroup */
class GriddedDataTimeSeries2PotentialCoefficients
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GriddedDataTimeSeries2PotentialCoefficients, PARALLEL, "Convert a griddedDataTimeSeries to a sequence of PotentialCoefficients files", Grid, TimeSeries)

/***********************************************/

void GriddedDataTimeSeries2PotentialCoefficients::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    std::vector<FileName> fileNameOut;
    std::string           nameTime, nameIndex, nameCount;
    FileName              fileNameIn;
    TimeSeriesPtr         timeSeries;
    std::vector<ExpressionVariablePtr> exprValue;
    ExpressionVariablePtr exprArea;
    KernelPtr             kernel;
    Double                GM, R;
    UInt                  minDegree, maxDegree;
    Bool                  useLeastSquares;

    readConfig(config, "outputfilesPotentialCoefficients", fileNameOut,     Config::MUSTSET,  "coeff_{loopTime:%D}.gfc", "for each epoch");
    readConfig(config, "variableLoopTime",                 nameTime,        Config::OPTIONAL, "loopTime", "variable with time of each epoch");
    readConfig(config, "variableLoopIndex",                nameIndex,       Config::OPTIONAL, "",         "variable with index of current epoch (starts with zero)");
    readConfig(config, "variableLoopCount",                nameCount,       Config::OPTIONAL, "",         "variable with total number of epochs");
    readConfig(config, "inputfileGriddedDataTimeSeries",   fileNameIn,      Config::MUSTSET,  "",         "");
    readConfig(config, "timeSeries",                       timeSeries,      Config::OPTIONAL, "",         "otherwise times from inputfile are used");
    readConfig(config, "value",                            exprValue,       Config::MUSTSET,  "data0",    "expression (variables: longitude, latitude, height, area, data0, data1, ...)");
    readConfig(config, "weight",                           exprArea,        Config::MUSTSET,  "area",     "expression to compute values (input columns are named data0, data1, ...)");
    readConfig(config, "kernel",                           kernel,          Config::MUSTSET,  "",         "kernel in which the grid values are given");
    readConfig(config, "minDegree",                        minDegree,       Config::DEFAULT,  "0",        "");
    readConfig(config, "maxDegree",                        maxDegree,       Config::MUSTSET,  "",         "");
    readConfig(config, "GM",                               GM,              Config::DEFAULT,  STRING_DEFAULT_GM, "Geocentric gravitational constant");
    readConfig(config, "R",                                R,               Config::DEFAULT,  STRING_DEFAULT_R,  "reference radius for potential coefficients");
    readConfig(config, "leastSquares",                     useLeastSquares, Config::DEFAULT,  "0",     "false: quadrature formular, true: least squares adjustment order by order");
    if(isCreateSchema(config)) return;

    logStatus<<"read gridded data time series <"<<fileNameIn<<">"<<Log::endl;
    InFileGriddedDataTimeSeries file(fileNameIn);
    GriddedData grid = file.grid();
    if(!grid.areas.size())
      grid.computeArea();
    MiscGriddedData::printStatistics(grid);
    std::vector<Time> times = file.times();
    if(timeSeries)
      times = timeSeries->times();

    // evaluate expressions
    // --------------------
    logStatus<<"calculate gridded data"<<Log::endl;
    VariableList varList;
    addDataVariables(grid, varList);
    exprArea ->simplify(varList);
    std::vector<Double> areas(grid.points.size());
    for(UInt i=0; i<grid.points.size(); i++)
    {
      evaluateDataVariables(grid, i, varList);
      areas.at(i) = exprArea->evaluate(varList);
    }

    // evaluate data at each epoch
    std::vector<std::vector<Double>> values(exprValue.size()*times.size(), std::vector<Double>(grid.points.size()));
    Single::forEach(times.size(), [&](UInt idEpoch)
    {
      Matrix data = file.data(times.at(idEpoch));
      grid.values.resize(data.columns());
      for(UInt k=0; k<data.columns(); k++)
        grid.values.at(k) = Vector(data.column(k));
      VariableList varList;
      addDataVariables(grid, varList);
      std::for_each(exprValue.begin(), exprValue.end(), [&](auto expr) {expr->simplify(varList);});
      for(UInt i=0; i<grid.points.size(); i++)
      {
        evaluateDataVariables(grid, i, varList);
        for(UInt k=0; k<exprValue.size(); k++)
          values.at(idEpoch*exprValue.size()+k).at(i) = exprValue.at(k)->evaluate(varList);
      }
    });
    grid.areas  = std::move(areas);
    grid.values = std::move(values);

    // spherical harmonic analysis
    // ---------------------------
    std::vector<SphericalHarmonics> harmonics = MiscGriddedData::analysisSphericalHarmonics(grid, kernel, minDegree, maxDegree, GM, R, useLeastSquares, comm);

    // write potential coefficients
    // ----------------------------
    if(Parallel::isMaster(comm))
    {
      logStatus<<"write "<<harmonics.size()<<" potential coefficient files"<<Log::endl;
      VariableList varList;
      if(!nameTime.empty())  varList.undefineVariable(nameTime);
      if(!nameIndex.empty()) varList.undefineVariable(nameIndex);
      if(!nameCount.empty()) varList.setVariable(nameCount, times.size());
      for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
      {
        if(!nameTime.empty())  varList.setVariable(nameTime, times.at(idEpoch).mjd());
        if(!nameIndex.empty()) varList.setVariable(nameIndex, idEpoch);
        for(UInt k=0; k<exprValue.size(); k++)
          writeFileSphericalHarmonics(fileNameOut.at(k)(varList), harmonics.at(idEpoch*exprValue.size()+k));
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
