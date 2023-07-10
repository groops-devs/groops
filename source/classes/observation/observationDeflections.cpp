/***********************************************/
/**
* @file observationDeflections.cpp
*
* @brief point measurements of xi and eta.
*
*
* @author Christian Pock
* @date 2012-05-30
*
*/
/***********************************************/

#include "base/import.h"
#include "base/planets.h"
#include "parser/dataVariables.h"
#include "files/fileGriddedData.h"
#include "classes/parametrizationGravity/parametrizationGravity.h"
#include "classes/observation/observation.h"
#include "classes/observation/observationDeflections.h"

/***********************************************/

ObservationDeflections::ObservationDeflections(Config &config)
{
  try
  {
    FileName              fileNameGrid;
    ExpressionVariablePtr exprXi, exprEta;
    ExpressionVariablePtr exprSigmaXi, exprSigmaEta;
    Double                a, f;

    renameDeprecatedConfig(config, "representation", "parametrizationGravity", date2time(2020, 6, 3));

    if(readConfigSequence(config, "rightHandSide", Config::MUSTSET, "", "input for observation vectors"))
    {
      readConfig(config, "inputfileGriddedData", fileNameGrid,   Config::MUSTSET,  "",                "");
      readConfig(config, "observationXi",        exprXi,         Config::MUSTSET,  "data0",           "North-South Deflections of the Vertical [rad]");
      readConfig(config, "observationEta",       exprEta,        Config::MUSTSET,  "data1",           "East-West Deflections of the Vertical  [rad]");
      readConfig(config, "sigmaXi",              exprSigmaXi,    Config::OPTIONAL, "sqrt(4*PI/area)", "accuracy, 1/sigma used as weighting");
      readConfig(config, "sigmaEta",             exprSigmaEta,   Config::OPTIONAL, "sqrt(4*PI/area)", "accuracy, 1/sigma used as weighting");
      readConfig(config, "referencefield",       referencefield, Config::OPTIONAL, "",                "");
      endSequence(config);
    }
    readConfig(config, "parametrizationGravity", parametrization, Config::MUSTSET,  "",   "");
    readConfig(config, "time",                   time,            Config::OPTIONAL, "",   "for reference field and parametrization");
    readConfig(config, "R",                      a,               Config::DEFAULT,  STRING_DEFAULT_GRS80_a, "reference radius for ellipsoid");
    readConfig(config, "inverseFlattening",      f,               Config::DEFAULT,  STRING_DEFAULT_GRS80_f, "reference flattening for ellipsoid, 0: sphere");
    readConfig(config, "blockingSize",           obsPerArc,       Config::DEFAULT, "100", "segementation of the obervations if designmatrix can't be build at once");
    if(isCreateSchema(config)) return;

    ellipsoid = Ellipsoid(a, f);

    GriddedData grid;
    readFileGriddedData(fileNameGrid, grid);
    points = grid.points;

    // evaluate expression
    // -------------------
    VariableList varList;
    addDataVariables(grid, varList);
    std::vector<ExpressionVariablePtr> expr({exprXi, exprEta, exprSigmaXi, exprSigmaEta});
    std::for_each(expr.begin(),  expr.end(),  [&](auto expr) {if(expr) expr->simplify(varList);});

    xi.resize(grid.points.size(), 0.);
    eta.resize(grid.points.size(), 0.);
    sigmasXi.resize(grid.points.size(), 1.);
    sigmasEta.resize(grid.points.size(), 1.);
    for(UInt i=0; i<grid.values.size(); i++)
    {
      evaluateDataVariables(grid, i, varList);
      if(exprXi)       xi.at(i)        = exprXi      ->evaluate(varList);
      if(exprEta)      eta.at(i)       = exprEta     ->evaluate(varList);
      if(exprSigmaXi)  sigmasXi.at(i)  = exprSigmaXi ->evaluate(varList);
      if(exprSigmaEta) sigmasEta.at(i) = exprSigmaEta->evaluate(varList);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ObservationDeflections::observation(UInt arcNo, Matrix &l, Matrix &A, Matrix &B)
{
  try
  {
    const UInt obsCount = observationCount(arcNo);

    l = Vector(2*obsCount);
    A = Matrix(2*obsCount, parameterCount());
    B = Matrix();
    for(UInt obsNo=0; obsNo<obsCount; obsNo++)
    {
      const Vector3d    pos     = points.at(idx(arcNo, obsNo));
      const Transform3d trf2neu = inverse(localNorthEastUp(pos, ellipsoid));
      Vector3d g;
      if(referencefield)
        g = trf2neu.transform(normalize(referencefield->gravity(time, pos)));

      l(2*obsNo+0, 0) = xi.at (idx(arcNo, obsNo)) - g.x();
      l(2*obsNo+1, 0) = eta.at(idx(arcNo, obsNo)) - g.y();

      Matrix G(3, parameterCount());
      parametrization->gravity(time, pos, G);
      matMult(-1./Planets::normalGravity(pos), trf2neu.matrix().row(0,2), G, A.row(2*obsNo, 2));

      l.row(2*obsNo+0) *= 1/sigmasXi.at (idx(arcNo, obsNo));
      l.row(2*obsNo+1) *= 1/sigmasEta.at(idx(arcNo, obsNo));
      A.row(2*obsNo+0) *= 1/sigmasXi.at (idx(arcNo, obsNo));
      A.row(2*obsNo+1) *= 1/sigmasEta.at(idx(arcNo, obsNo));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
