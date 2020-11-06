/***********************************************/
/**
* @file observationTerrestrial.cpp
*
* @brief Terrestrial observations (point measurements).
* Observed function values of Gravityfield.
*
* @author Torsten Mayer-Guerr
* @date 2005-03-22
*
*/
/***********************************************/

#include "base/import.h"
#include "parser/dataVariables.h"
#include "files/fileGriddedData.h"
#include "classes/kernel/kernel.h"
#include "classes/gravityfield/gravityfield.h"
#include "classes/parametrizationGravity/parametrizationGravity.h"
#include "classes/observation/observation.h"
#include "classes/observation/observationTerrestrial.h"

/***********************************************/

ObservationTerrestrial::ObservationTerrestrial(Config &config)
{
  try
  {
    FileName              fileNameGrid;
    ExpressionVariablePtr exprValue, exprSigma;

    renameDeprecatedConfig(config, "representation", "parametrizationGravity", date2time(2020, 6, 3));

    if(readConfigSequence(config, "rightHandSide", Config::MUSTSET, "", "input for observation vectors"))
    {
      readConfig(config, "inputfileGriddedData", fileNameGrid,   Config::MUSTSET,  "",                "");
      readConfig(config, "observation",          exprValue,      Config::MUSTSET,  "data0",           "[SI units]");
      readConfig(config, "sigma",                exprSigma,      Config::OPTIONAL, "sqrt(4*PI/area)", "accuracy, 1/sigma used as weighting");
      readConfig(config, "referencefield",       referencefield, Config::DEFAULT,  "",                "");
      endSequence(config);
    }
    readConfig(config, "kernel",                 kernel,          Config::MUSTSET,  "",    "type of observations");
    readConfig(config, "parametrizationGravity", parametrization, Config::MUSTSET,  "",    "");
    readConfig(config, "time",                   time,            Config::OPTIONAL, "",    "for reference field and parametrization");
    readConfig(config, "blockingSize",           obsPerArc,       Config::DEFAULT,  "100", "segementation of the obervations if designmatrix can't be build at once");
    if(isCreateSchema(config)) return;

    GriddedData grid;
    readFileGriddedData(fileNameGrid, grid);
    points = grid.points;

    // evaluate expression
    // -------------------
    auto varList = config.getVarList();
    std::set<std::string> usedVariables;
    if(exprValue) exprValue->usedVariables(varList, usedVariables);
    if(exprSigma) exprSigma->usedVariables(varList, usedVariables);
    addDataVariables(grid, varList, usedVariables);
    if(exprValue) exprValue->simplify(varList);
    if(exprSigma) exprSigma->simplify(varList);

    values.resize(grid.points.size(), 0.);
    sigmas.resize(grid.points.size(), 1.);
    for(UInt i=0; i<values.size(); i++)
    {
      evaluateDataVariables(grid, i, varList);
      if(exprValue) values.at(i) = exprValue->evaluate(varList);
      if(exprSigma) sigmas.at(i) = exprSigma->evaluate(varList);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ObservationTerrestrial::observation(UInt arcNo, Matrix &l, Matrix &A, Matrix &B)
{
  try
  {
    const UInt obsCount = observationCount(arcNo);

    l = Vector(obsCount);
    A = Matrix(obsCount, parameterCount());
    B = Matrix();
    for(UInt obsNo=0; obsNo<obsCount; obsNo++)
    {
      l(obsNo, 0) = values.at(idx(arcNo, obsNo)) - referencefield->field(time, points.at(idx(arcNo, obsNo)), *kernel);
      parametrization->field(time, points.at(idx(arcNo, obsNo)), *kernel, A.row(obsNo));
      if(sigmas.at(idx(arcNo, obsNo)) != 1.)
      {
        l.row(obsNo) *= 1/sigmas.at(idx(arcNo, obsNo));
        A.row(obsNo) *= 1/sigmas.at(idx(arcNo, obsNo));
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
