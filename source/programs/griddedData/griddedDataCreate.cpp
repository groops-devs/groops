/***********************************************/
/**
* @file griddedDataCreate.cpp
*
* @brief Create a grid an write it to file.
*
* @author Torsten Mayer-Guerr
* @date 2001-08-06
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program creates a \configClass{grid}{gridType} and writes it to \configFile{outputfileGrid}{griddedData}.
The grid is expressed as ellipsoidal coordinates (longitude, latitude, height)
based on a reference ellipsoid with parameters \config{R} and \config{inverseFlattening}.
Extra \config{value} columns can be appended using expressions
with the common \reference{data variables}{general.parser:dataVariables} for gridded data.
)";

/***********************************************/

#include "programs/program.h"
#include "parser/dataVariables.h"
#include "files/fileGriddedData.h"
#include "classes/grid/grid.h"
#include "misc/miscGriddedData.h"

/***** CLASS ***********************************/

/** @brief Create a grid an write it to file.
* @ingroup programsGroup */
class GriddedDataCreate
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GriddedDataCreate, SINGLEPROCESS, "create a grid an write it to file", Grid)
GROOPS_RENAMED_PROGRAM(Grid2File, GriddedDataCreate, date2time(2020, 01, 28))

/***********************************************/

void GriddedDataCreate::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameGrid;
    GridPtr  gridPtr;
    Double   a, f;
    std::vector<ExpressionVariablePtr> exprValues;

    readConfig(config, "outputfileGrid",    fileNameGrid, Config::MUSTSET,  "",    "");
    readConfig(config, "grid",              gridPtr,      Config::MUSTSET,  "",    "");
    readConfig(config, "R",                 a,            Config::DEFAULT,  STRING_DEFAULT_GRS80_a, "reference radius for ellipsoidal coordinates on output");
    readConfig(config, "inverseFlattening", f,            Config::DEFAULT,  STRING_DEFAULT_GRS80_f, "reference flattening for ellipsoidal coordinates on output, 0: spherical coordinates");
    readConfig(config, "value",             exprValues,   Config::OPTIONAL, "0.0", "expression (variables as 'longitude', 'height', 'area' are taken from the gridded data)");
    if(isCreateSchema(config)) return;

    logStatus<<"create grid"<<Log::endl;
    GriddedData grid(Ellipsoid(a,f), gridPtr->points(), gridPtr->areas(), {});

    if(exprValues.size())
    {
      VariableList varList;
      addDataVariables(grid, varList);
      for(auto exprValue : exprValues)
        exprValue->simplify(varList);

      grid.values.resize(exprValues.size(), std::vector<Double>(grid.points.size(), 0));
      for(UInt i=0; i<grid.points.size(); i++)
      {
        evaluateDataVariables(grid, i, varList);
        for(UInt k=0; k<exprValues.size(); k++)
          grid.values.at(k).at(i) = exprValues.at(k)->evaluate(varList);
      }
    }

    logStatus<<"writing grid to file <"<<fileNameGrid<<">"<<Log::endl;
    writeFileGriddedData(fileNameGrid, grid);
    MiscGriddedData::printStatistics(grid);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
