/***********************************************/
/**
* @file matrix2GriddedData.cpp
*
* @brief Read gridded data from matrix.
*
* @author Torsten Mayer-Guerr
* @date 2006-03-17
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program reads a \file{matrix file}{matrix} with data in columns
and convert into \file{gridded data}{griddedData}.
The input columns are enumerated by \verb|data0|,~\verb|data1|,~\ldots,
see~\reference{dataVariables}{general.parser:dataVariables}.
)";

/***********************************************/

#include "programs/program.h"
#include "parser/dataVariables.h"
#include "files/fileMatrix.h"
#include "files/fileGriddedData.h"
#include "misc/miscGriddedData.h"

/***** CLASS ***********************************/

/** @brief Read gridded data from matrix.
* @ingroup programsGroup */
class Matrix2GriddedData
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(Matrix2GriddedData, SINGLEPROCESS, "Read gridded data from matrix", Grid, Matrix)
GROOPS_RENAMED_PROGRAM(GridMatrix2GriddedData, Matrix2GriddedData, date2time(2020, 02, 12))

/***********************************************/

void Matrix2GriddedData::run(Config &config)
{
  try
  {
    FileName              fileNameOut, fileNameIn;
    ExpressionVariablePtr exprLon, exprLat, exprHeight;
    ExpressionVariablePtr exprX, exprY, exprZ;
    ExpressionVariablePtr exprArea;
    std::vector<ExpressionVariablePtr> exprValues;
    Bool                  sortPoints, computeArea;
    Double                a, f;
    std::string           choice;

    readConfig(config, "outputfileGriddedData", fileNameOut, Config::MUSTSET,  "", "");
    readConfig(config, "inputfileMatrix",       fileNameIn,  Config::MUSTSET,  "", "");
    if(readConfigChoice(config, "points", choice,  Config::MUSTSET, "", ""))
    {
      if(readConfigChoiceElement(config, "ellipsoidal", choice, ""))
      {
        readConfig(config, "longitude", exprLon,    Config::MUSTSET,  "data0", "expression");
        readConfig(config, "latitude",  exprLat,    Config::MUSTSET,  "data1", "expression");
        readConfig(config, "height",    exprHeight, Config::MUSTSET,  "data2", "expression");
      }
      if(readConfigChoiceElement(config, "cartesian", choice, ""))
      {
        readConfig(config, "x", exprX, Config::MUSTSET,  "data0", "expression");
        readConfig(config, "y", exprY, Config::MUSTSET,  "data1", "expression");
        readConfig(config, "z", exprZ, Config::MUSTSET,  "data2", "expression");
      }
      endChoice(config);
    }
    readConfig(config, "area",              exprArea,    Config::OPTIONAL, "", "expression (e.g. deltaL*2*sin(deltaB/2)*cos(data1/RHO))");
    readConfig(config, "value",             exprValues,  Config::OPTIONAL, "data3", "expression");
    readConfig(config, "sortPoints",        sortPoints,  Config::DEFAULT,  "0", "sort from north/west to south east");
    readConfig(config, "computeArea",       computeArea, Config::DEFAULT,  "0", "the area can be computed automatically for rectangular grids");
    readConfig(config, "R",                 a,           Config::DEFAULT,  STRING_DEFAULT_GRS80_a, "reference radius for ellipsoidal coordinates");
    readConfig(config, "inverseFlattening", f,           Config::DEFAULT,  STRING_DEFAULT_GRS80_f, "reference flattening for ellipsoidal coordinates");
    if(isCreateSchema(config)) return;

    // =====================================================

    logStatus<<"Read matrix from file <"<<fileNameIn<<">"<<Log::endl;
    Matrix A;
    readFileMatrix(fileNameIn, A);
    logInfo<<"  matrix: "<<A.rows()<<" x "<<A.columns()<<Log::endl;

    std::vector<ExpressionVariablePtr> expr({exprX, exprY, exprZ, exprArea}); // collect all expr
    expr.insert(expr.end(), exprValues.begin(), exprValues.end());            // collect all expr
    auto varList = config.getVarList();
    std::set<std::string> usedVariables;
    std::for_each(expr.begin(), expr.end(), [&](auto expr) {if(expr) expr->usedVariables(varList, usedVariables);});
    addDataVariables(A, varList, usedVariables);
    std::for_each(exprValues.begin(), exprValues.end(), [&](auto expr) {if(expr) expr->simplify(varList);});

    // create grid
    // -----------
    logStatus<<"create grid"<<Log::endl;
    GriddedData grid;
    grid.ellipsoid = Ellipsoid(a,f);
    grid.values.resize(exprValues.size());
    for(UInt i=0; i<A.rows(); i++)
    {
      evaluateDataVariables(A, i, varList);
      if(exprLon)
        grid.points.push_back(grid.ellipsoid(Angle(DEG2RAD*exprLon->evaluate(varList)), Angle(DEG2RAD*exprLat->evaluate(varList)), exprHeight->evaluate(varList)));
      else
        grid.points.push_back(Vector3d(exprX->evaluate(varList), exprY->evaluate(varList), exprZ->evaluate(varList)));
      if(exprArea)
        grid.areas.push_back( exprArea->evaluate(varList) );
      for(UInt k=0; k<exprValues.size(); k++)
        grid.values.at(k).push_back( exprValues.at(k)->evaluate(varList) );
    } // for(i)

    // =====================================================

    if(sortPoints)
    {
      logStatus<<"sort points"<<Log::endl;
      grid.sort();
    }

    if(computeArea)
    {
      logStatus<<"compute area"<<Log::endl;
      if(!grid.computeArea())
        logWarning<<"Compute areas only possible with regular rectangular grids"<<Log::endl;
    }

    // save grid
    // ---------
    logStatus<<"save grid <"<<fileNameOut<<">"<<Log::endl;
    writeFileGriddedData(fileNameOut, grid);
    MiscGriddedData::printStatistics(grid);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
