/***********************************************/
/**
* @file observationStationLoading.cpp
*
* @brief Loading from station observations.
*
* @author Torsten Mayer-Guerr
* @date 2014-07-17
*
*/
/***********************************************/

#include "base/import.h"
#include "base/planets.h"
#include "parser/dataVariables.h"
#include "files/fileMatrix.h"
#include "files/fileGriddedData.h"
#include "classes/gravityfield/gravityfield.h"
#include "classes/parametrizationGravity/parametrizationGravity.h"
#include "classes/observation/observation.h"
#include "classes/observation/observationStationLoading.h"

/***********************************************/

ObservationStationLoading::ObservationStationLoading(Config &config)
{
  try
  {
    FileName              fileNameGrid;
    ExpressionVariablePtr exprNorth, exprEast, exprUp;
    ExpressionVariablePtr exprSigmaNorth, exprSigmaEast, exprSigmaUp;
    FileName              deformationName, potentialName;
    Double                a, f;

    renameDeprecatedConfig(config, "representation", "parametrizationGravity", date2time(2020, 6, 3));

    if(readConfigSequence(config, "rightHandSide", Config::MUSTSET, "", "input for observation vectors"))
    {
      readConfig(config, "inputfileGriddedData", fileNameGrid,   Config::MUSTSET,  "",      "station positions with displacement data");
      readConfig(config, "observationNorth",     exprNorth,      Config::MUSTSET,  "data0", "displacement [m]");
      readConfig(config, "observationEast",      exprEast,       Config::MUSTSET,  "data1", "displacement [m]");
      readConfig(config, "observationUp",        exprUp,         Config::MUSTSET,  "data2", "displacement [m]");
      readConfig(config, "sigmaNorth",           exprSigmaNorth, Config::OPTIONAL, "data3", "accuracy, 1/sigma used as weighting");
      readConfig(config, "sigmaEast",            exprSigmaEast,  Config::OPTIONAL, "data4", "accuracy, 1/sigma used as weighting");
      readConfig(config, "sigmaUp",              exprSigmaUp,    Config::OPTIONAL, "data5", "accuracy, 1/sigma used as weighting");
      readConfig(config, "inGlobalFrame",        inGlobalFrame,  Config::DEFAULT,  "0",     "obs/sigmas given in global x,y,z frame instead of north,east,up");
      readConfig(config, "referencefield",       referencefield, Config::OPTIONAL, "",      "");
      endSequence(config);
    }
    readConfig(config, "time",                   time,                Config::OPTIONAL, "",  "for reference field and parametrization");
    readConfig(config, "parametrizationGravity", parametrization,     Config::MUSTSET,  "",  "of loading (+defo) potential");
    readConfig(config, "estimateTranslation",    estimateTranslation, Config::DEFAULT,  "1", "coordinate center");
    readConfig(config, "estimateScale",          estimateScale,       Config::DEFAULT,  "1", "scale factor of position");
    readConfig(config, "estimateRotation",       estimateRotation,    Config::DEFAULT,  "1", "rotation");
    readConfig(config, "inputfileDeformationLoadLoveNumber", deformationName, Config::MUSTSET,  "{groopsDataDir}/loading/deformationLoveNumbers_CM_Gegout97.txt", "");
    readConfig(config, "inputfilePotentialLoadLoveNumber",   potentialName,   Config::OPTIONAL, "{groopsDataDir}/loading/loadLoveNumbers_Gegout97.txt", "if full potential is given and not only loading potential");
    readConfig(config, "R",                      a,                   Config::DEFAULT,  STRING_DEFAULT_GRS80_a, "reference radius for ellipsoid");
    readConfig(config, "inverseFlattening",      f,                   Config::DEFAULT,  STRING_DEFAULT_GRS80_f, "reference flattening for ellipsoid, 0: sphere");
    if(isCreateSchema(config)) return;

    ellipsoid = Ellipsoid(a, f);

    GriddedData grid;
    readFileGriddedData(fileNameGrid, grid);
    points = grid.points;

    // evaluate expression
    // -------------------
    VariableList varList;
    addDataVariables(grid, varList);
    std::vector<ExpressionVariablePtr> expr({exprNorth, exprEast, exprUp, exprSigmaNorth, exprSigmaEast, exprSigmaUp});
    std::for_each(expr.begin(), expr.end(), [&](auto expr) {if(expr) expr->simplify(varList);});

    valuesNorth.resize(grid.points.size(), 0.);
    valuesEast.resize (grid.points.size(), 0.);
    valuesUp.resize   (grid.points.size(), 0.);
    sigmasNorth.resize(grid.points.size(), 1.);
    sigmasEast.resize (grid.points.size(), 1.);
    sigmasUp.resize   (grid.points.size(), 1.);
    for(UInt i=0; i<grid.values.size(); i++)
    {
      evaluateDataVariables(grid, i, varList);
      if(exprNorth)      valuesNorth.at(i) = exprNorth     ->evaluate(varList);
      if(exprEast)       valuesEast.at(i)  = exprEast      ->evaluate(varList);
      if(exprUp)         valuesUp.at(i)    = exprUp        ->evaluate(varList);
      if(exprSigmaNorth) sigmasNorth.at(i) = exprSigmaNorth->evaluate(varList);
      if(exprSigmaEast)  sigmasEast.at(i)  = exprSigmaEast ->evaluate(varList);
      if(exprSigmaUp)    sigmasUp.at(i)    = exprSigmaUp   ->evaluate(varList);
    }

    // load love numbers
    // -----------------
    Matrix love;
    readFileMatrix(deformationName, love);
    hn = love.column(0);
    ln = love.column(1);

    // models contain the total mass (loading mass & deformation mass effect)
    if(!potentialName.empty())
    {
      Vector kn;
      readFileMatrix(potentialName, kn);
      for(UInt n=2; n<std::min(kn.rows(), hn.rows()); n++)
        hn(n) /= (1.+kn(n));
      for(UInt n=2; n<std::min(kn.rows(), ln.rows()); n++)
        ln(n) /= (1.+kn(n));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ObservationStationLoading::parameterName(std::vector<ParameterName> &name) const
{
  try
  {
    parametrization->parameterName(name);
    if(estimateTranslation)
    {
      name.push_back(ParameterName("", "translation.x"));
      name.push_back(ParameterName("", "translation.y"));
      name.push_back(ParameterName("", "translation.z"));
    }
    if(estimateScale)
    {
      name.push_back(ParameterName("", "scale"));
    }
    if(estimateRotation)
    {
      name.push_back(ParameterName("", "rotation.x"));
      name.push_back(ParameterName("", "rotation.y"));
      name.push_back(ParameterName("", "rotation.z"));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ObservationStationLoading::observation(UInt /*arcNo*/, Matrix &l, Matrix &A, Matrix &B)
{
  try
  {
    // Local (north, east, down) -> global TRF
    // -------------------------------------
    std::vector<Transform3d> trf2neu(points.size());
      for(UInt i=0; i<points.size(); i++)
        trf2neu.at(i) = inverse(localNorthEastUp(points.at(i), ellipsoid));

    // normal gravity
    // --------------
    std::vector<Double> gravity(points.size());
    for(UInt i=0; i<points.size(); i++)
      gravity.at(i) = Planets::normalGravity(points.at(i));

    // observations
    // ------------
    l = Matrix(3*points.size(), 1);
    for(UInt i=0; i<points.size(); i++)
    {
      const Vector3d obs = Vector3d(valuesNorth.at(i), valuesEast.at(i), valuesUp.at(i));
      l(3*i+0, 0) = obs.x();
      l(3*i+1, 0) = obs.y();
      l(3*i+2, 0) = obs.z();
    }

    // reduce reference displacement
    // -----------------------------
    if(referencefield)
    {
      std::vector<std::vector<Vector3d>> referenceDefo(points.size(), std::vector<Vector3d>(1));
      referencefield->deformation({time}, points, gravity, hn, ln, referenceDefo);
      for(UInt i=0; i<points.size(); i++)
      {
        Vector3d d = referenceDefo.at(i).at(0);
        if(!inGlobalFrame)
          d = trf2neu.at(i).transform(d);
        l(3*i+0, 0) -= d.x();
        l(3*i+1, 0) -= d.y();
        l(3*i+2, 0) -= d.z();
      }
    }

    // design matrix A
    // ---------------
    A = Matrix(3*points.size(), parameterCount());
    B = Matrix();
    for(UInt i=0; i<points.size(); i++)
    {
      parametrization->deformation(time, points.at(i), gravity.at(i), hn, ln, A.slice(3*i, 0, 3, parametrization->parameterCount()));

      // helmert transformation
      UInt idx = parametrization->parameterCount();
      if(estimateTranslation)
      {
        A(3*i+0, idx+0) = 1.;
        A(3*i+1, idx+1) = 1.;
        A(3*i+2, idx+2) = 1.;
        idx += 3;
      }
      if(estimateScale)
      {
        A(3*i+0, idx+3) =  points.at(i).x();
        A(3*i+1, idx+3) =  points.at(i).y();
        A(3*i+2, idx+3) =  points.at(i).z();
        idx++;
      }
      if(estimateRotation)
      {
        A(3*i+0, idx+5) =  points.at(i).z();
        A(3*i+0, idx+6) = -points.at(i).y();
        A(3*i+1, idx+4) = -points.at(i).z();
        A(3*i+1, idx+6) =  points.at(i).x();
        A(3*i+2, idx+4) =  points.at(i).y();
        A(3*i+2, idx+5) = -points.at(i).x();
        idx += 3;
      }

      // rotate into local system
      if(!inGlobalFrame)
        copy(trf2neu.at(i).matrix()*A.row(3*i,3), A.row(3*i,3));
    }

    // Decorrelation
    // -------------
    for(UInt i=0; i<points.size(); i++)
    {
      l.row(3*i+0) *= 1/sigmasNorth.at(i);
      l.row(3*i+1) *= 1/sigmasEast.at(i);
      l.row(3*i+2) *= 1/sigmasUp.at(i);
      A.row(3*i+0) *= 1/sigmasNorth.at(i);
      A.row(3*i+1) *= 1/sigmasEast.at(i);
      A.row(3*i+2) *= 1/sigmasUp.at(i);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
