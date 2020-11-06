/***********************************************/
/**
* @file dataVariables.cpp
*
* @brief Create variables to evaluate a matrix/grid/timeSeries.
*
* @author Torsten Mayer-Guerr
* @date 2018-06-18
*
*/
/***********************************************/

#include "base/importStd.h"
#include "base/time.h"
#include "base/griddedData.h"
#include "parser/stringParser.h"
#include "parser/expressionParser.h"
#include "parser/dataVariables.h"

/***********************************************/

void addTimeVariables(VariableList &varList)
{
  try
  {
    addVariable("loopTime",      varList);
    addVariable("loopTimeStart", varList);
    addVariable("loopTimeEnd",   varList);
    addVariable("index",         varList);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void evaluateTimeVariables(UInt index, const Time &timeStart, const Time &timeEnd, VariableList &varList)
{
  try
  {
    varList["loopTime"]->setValue(timeStart.mjd());
    varList["loopTimeStart"]->setValue(timeStart.mjd());
    varList["loopTimeEnd"]->setValue(timeEnd.mjd());
    varList["index"]->setValue(static_cast<Double>(index));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void addDataVariables(const std::string &prefix, const_MatrixSliceRef data, VariableList &varList, const std::set<std::string> &usedName)
{
  try
  {
    auto used = [&](const std::string &name) {return (usedName.find(name) != usedName.end());};

    Double median = NAN_EXPR;
    Double mad    = NAN_EXPR; // Median Absolute Deviation
    if(data.rows() && (used(prefix+"median") || used(prefix+"mad")))
    {
      median = ::median(data);
      if(used(prefix+"mad"))
      {
        Matrix tmp = data;
        for(UInt i=0; i<tmp.rows(); i++)
          for(UInt k=0; k<tmp.columns(); k++)
            tmp(i,k) = std::abs(tmp(i,k)-median);
        mad = ::median(tmp);
      }
    }

    Double step = NAN_EXPR;
    if((data.size()>1) && used(prefix+"step"))
    {
      std::vector<Double> tmp = flatten(data);
      std::adjacent_difference(tmp.begin(), tmp.end(), tmp.begin());
      std::transform(tmp.begin(), tmp.end(), tmp.begin(), [](const Double &x) {return std::fabs(x);});
      step = *std::min_element(tmp.begin(), tmp.end());
    }

    addVariable(prefix, varList);
    if(used(prefix+"count"))  addVariable(prefix+"count",  static_cast<Double>(data.size()), varList);
    if(used(prefix+"min"))    addVariable(prefix+"min",    min(data),                        varList);
    if(used(prefix+"max"))    addVariable(prefix+"max",    max(data),                        varList);
    if(used(prefix+"mean"))   addVariable(prefix+"mean",   mean(data),                       varList);
    if(used(prefix+"rms"))    addVariable(prefix+"rms",    rootMeanSquare(data),             varList);
    if(used(prefix+"std"))    addVariable(prefix+"std",    standardDeviation(data),          varList);
    if(used(prefix+"median")) addVariable(prefix+"median", median,                           varList);
    if(used(prefix+"mad"))    addVariable(prefix+"mad",    mad,                              varList);
    if(used(prefix+"step"))   addVariable(prefix+"step",   step,                             varList);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void addDataVariables(const std::string &prefix, const std::vector<Time> &times, VariableList &varList, const std::set<std::string> &usedName)
{
  try
  {
    Vector data(times.size());
    for(UInt i=0; i<times.size(); i++)
      data(i) = times.at(i).mjd();
    addDataVariables(prefix, data, varList, usedName);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void addDataVariables(const std::string &prefix, const_MatrixSliceRef data, const_MatrixSliceRef weight, VariableList &varList, const std::set<std::string> &usedName)
{
  try
  {
    addDataVariables(prefix, data, varList, usedName);
    if(!weight.size())
      return;

    auto used = [&](const std::string &name) {return (usedName.find(name) != usedName.end());};

    Double w     = sum(weight);
    Double wmean = inner(data, weight)/w;
    Double wrms  = 0;
    for(UInt i=0; i<data.rows(); i++)
      for(UInt k=0; k<data.columns(); k++)
        wrms += weight(i,k)/w * data(i,0) * data(i,k);

    if(used(prefix+"wmean")) addVariable(prefix+"wmean", wmean, varList);
    if(used(prefix+"wrms"))  addVariable(prefix+"wrms",  std::sqrt(wrms), varList);
    if(used(prefix+"wstd"))  addVariable(prefix+"wstd",  std::sqrt(data.size()/(data.size()-1.)*(wrms-wmean*wmean)), varList);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void addDataVariables(const_MatrixSliceRef data, VariableList &varList, const std::set<std::string> &usedName)
{
  try
  {
    addVariable("index", varList);
    for(UInt i=0; i<data.columns(); i++)
      addDataVariables("data"+i%"%i"s, data.column(i), varList, usedName);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void evaluateDataVariables(const_MatrixSliceRef data, UInt row, VariableList &varList)
{
  try
  {
    varList["index"]->setValue(static_cast<Double>(row)); // index
    for(UInt i=0; i<data.columns(); i++)
      varList["data"+i%"%i"s]->setValue(data(row,i)); // "data(i)"
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void undefineDataVariables(const_MatrixSliceRef data, VariableList &varList)
{
  try
  {
    varList["index"]->setUndefined(); // index
    for(UInt i=0; i<data.columns(); i++)
      varList["data"+i%"%i"s]->setUndefined(); // "data(i)"
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void addDataVariables(const GriddedData &grid, VariableList &varList, const std::set<std::string> &usedName)
{
  try
  {
    addVariable("longitude",  varList);
    addVariable("latitude",   varList);
    addVariable("height",     varList);
    addVariable("cartesianX", varList);
    addVariable("cartesianY", varList);
    addVariable("cartesianZ", varList);
    addVariable("area",       varList);
    addVariable("index",      varList);
    for(UInt i=0; i<grid.values.size(); i++)
      addDataVariables("data"+i%"%i"s, Vector(grid.values.at(i)), Vector(grid.areas), varList, usedName);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void evaluateDataVariables(const GriddedData &grid, UInt row, VariableList &varList)
{
  try
  {
    Angle  L, B;
    Double h;
    grid.ellipsoid(grid.points.at(row), L, B, h);
    varList["longitude"]->setValue(Double(L)*RAD2DEG);
    varList["latitude"]->setValue(Double(B)*RAD2DEG);
    varList["height"]->setValue(h);
    varList["cartesianX"]->setValue(grid.points.at(row).x());
    varList["cartesianY"]->setValue(grid.points.at(row).y());
    varList["cartesianZ"]->setValue(grid.points.at(row).z());
    varList["area"]->setValue((grid.areas.size() ? grid.areas.at(row) : 0));
    varList["index"]->setValue(static_cast<Double>(row));
    for(UInt i=0; i<grid.values.size(); i++)
      varList["data"+i%"%i"s]->setValue(grid.values.at(i).at(row));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void undefineDataVariables(const GriddedData &grid, VariableList &varList)
{
  try
  {
    varList["longitude"]->setUndefined();
    varList["latitude"]->setUndefined();
    varList["height"]->setUndefined();
    varList["cartesianX"]->setUndefined();
    varList["cartesianY"]->setUndefined();
    varList["cartesianZ"]->setUndefined();
    varList["area"]->setUndefined();
    varList["index"]->setUndefined();
    for(UInt i=0; i<grid.values.size(); i++)
      varList["data"+i%"%i"s]->setUndefined();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void addDataVariables(const GriddedDataRectangular &grid, VariableList &varList, const std::set<std::string> &usedName)
{
  try
  {
    auto used = [&](const std::string &name) {return (usedName.find(name) != usedName.end());};

    std::vector<Angle>  lambda, phi;
    std::vector<Double> radius;
    std::vector<Double> dLambda, dPhi;
    grid.geocentric(lambda, phi, radius, dLambda, dPhi);
    const UInt rows = phi.size();
    const UInt cols = lambda.size();

    // compute area
    for(UInt z=0; z<rows; z++)
      dPhi.at(z) *= std::cos(phi.at(z));
    const Double totalArea = std::accumulate(dLambda.begin(), dLambda.end(), Double(0.)) * std::accumulate(dPhi.begin(), dPhi.end(), Double(0.));

    addVariable("longitude",  varList);
    addVariable("latitude",   varList);
    addVariable("height",     varList);
    addVariable("cartesianX", varList);
    addVariable("cartesianY", varList);
    addVariable("cartesianZ", varList);
    addVariable("index",      varList);
    for(UInt i=0; i<grid.values.size(); i++)
    {
      const std::string prefix = "data"+i%"%i"s;
      addDataVariables(prefix, grid.values.at(i), varList, usedName);

      Double wmean = 0;
      Double wrms  = 0;
      for(UInt z=0; z<rows; z++)
        for(UInt s=0; s<cols; s++)
        {
          const Double w = dLambda.at(s)*dPhi.at(z)/totalArea;
          wmean += w * grid.values.at(i)(z,s);
          wrms  += w * grid.values.at(i)(z,s) * grid.values.at(i)(z,s);
        }

      if(used(prefix+"wmean")) addVariable(prefix+"wmean", wmean, varList);
      if(used(prefix+"wrms"))  addVariable(prefix+"wrms",  std::sqrt(wrms), varList);
      if(used(prefix+"wstd"))  addVariable(prefix+"wstd",  std::sqrt(rows*cols/(rows*cols-1.)*(wrms-wmean*wmean)), varList);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void evaluateDataVariables(const GriddedDataRectangular &grid, UInt row, UInt col, VariableList &varList)
{
  try
  {
    // set values of data variables
    const Vector3d p = grid.ellipsoid(grid.longitudes.at(col), grid.latitudes.at(row), grid.heights.at(row));
    varList["longitude"]->setValue(grid.longitudes.at(col)*RAD2DEG);
    varList["latitude"]->setValue(grid.latitudes.at(row)*RAD2DEG);
    varList["height"]->setValue(grid.heights.at(row));
    varList["cartesianX"]->setValue(p.x());
    varList["cartesianY"]->setValue(p.y());
    varList["cartesianZ"]->setValue(p.z());
    varList["index"]->setValue(static_cast<Double>(row*grid.longitudes.size()+col));
    for(UInt i=0; i<grid.values.size(); i++)
      varList["data"+i%"%i"s]->setValue(grid.values.at(i)(row,col)); // "data(i)"
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void undefineDataVariables(const GriddedDataRectangular &grid, VariableList &varList)
{
  try
  {
    varList["longitude"]->setUndefined();
    varList["latitude"]->setUndefined();
    varList["height"]->setUndefined();
    varList["cartesianX"]->setUndefined();
    varList["cartesianY"]->setUndefined();
    varList["cartesianZ"]->setUndefined();
    varList["index"]->setUndefined();
    for(UInt i=0; i<grid.values.size(); i++)
      varList["data"+i%"%i"s]->setUndefined(); // "data(i)"
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/
