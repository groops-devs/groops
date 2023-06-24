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
    varList.undefineVariable("loopTime");
    varList.undefineVariable("loopTimeStart");
    varList.undefineVariable("loopTimeEnd");
    varList.undefineVariable("index");
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
    varList.setVariable("loopTime", timeStart.mjd());
    varList.setVariable("loopTimeStart", timeStart.mjd());
    varList.setVariable("loopTimeEnd", timeEnd.mjd());
    varList.setVariable("index", static_cast<Double>(index));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

static Double step(const_MatrixSliceRef data)
{
  if(data.size() < 2)
    return NAN_EXPR;
  std::vector<Double> tmp = flatten(data);
  std::adjacent_difference(tmp.begin(), tmp.end(), tmp.begin());
  std::transform(tmp.begin(), tmp.end(), tmp.begin(), [](const Double &x) {return std::fabs(x);});
  return *std::min_element(tmp.begin(), tmp.end());
}

/***********************************************/

void addDataVariables(const std::string &prefix, const_MatrixSliceRef data, VariableList &varList)
{
  try
  {
    class Func : public ExpressionVariable::Func
    {
      MatrixSlice data;
      Double (*f)(const_MatrixSliceRef);
    public:
      Func(const_MatrixSliceRef data, Double (*f)(const_MatrixSliceRef)) : data(data), f(f) {}
      virtual ~Func() {}
      Double operator()() const override {return f(data);}
    };

    varList.undefineVariable(prefix);
    varList.setVariable(prefix+"count", static_cast<Double>(data.size()));
    varList.addVariable(std::make_shared<ExpressionVariable>(prefix+"min",    std::make_shared<Func>(data, min)));
    varList.addVariable(std::make_shared<ExpressionVariable>(prefix+"max",    std::make_shared<Func>(data, max)));
    varList.addVariable(std::make_shared<ExpressionVariable>(prefix+"sum",    std::make_shared<Func>(data, sum)));
    varList.addVariable(std::make_shared<ExpressionVariable>(prefix+"mean",   std::make_shared<Func>(data, mean)));
    varList.addVariable(std::make_shared<ExpressionVariable>(prefix+"rms",    std::make_shared<Func>(data, rootMeanSquare)));
    varList.addVariable(std::make_shared<ExpressionVariable>(prefix+"std",    std::make_shared<Func>(data, standardDeviation)));
    varList.addVariable(std::make_shared<ExpressionVariable>(prefix+"median", std::make_shared<Func>(data, median)));
    varList.addVariable(std::make_shared<ExpressionVariable>(prefix+"mad",    std::make_shared<Func>(data, medianAbsoluteDeviation)));
    varList.addVariable(std::make_shared<ExpressionVariable>(prefix+"step",   std::make_shared<Func>(data, step)));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void addDataVariables(const std::string &prefix, const std::vector<Time> &times, VariableList &varList)
{
  try
  {
    Vector data(times.size());
    for(UInt i=0; i<times.size(); i++)
      data(i) = times.at(i).mjd();
    addDataVariables(prefix, data, varList);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void addDataVariables(const std::string &prefix, const_MatrixSliceRef data, const_MatrixSliceRef weight, VariableList &varList)
{
  try
  {
    addDataVariables(prefix, data, varList);
    if(!weight.size())
      return;

    Double w     = sum(weight);
    Double wmean = inner(data, weight)/w;
    Double wrms  = 0;
    for(UInt i=0; i<data.rows(); i++)
      for(UInt k=0; k<data.columns(); k++)
        wrms += weight(i,k)/w * data(i,0) * data(i,k);

    varList.setVariable(prefix+"wmean", wmean);
    varList.setVariable(prefix+"wrms",  std::sqrt(wrms));
    varList.setVariable(prefix+"wstd",  std::sqrt(data.size()/(data.size()-1.)*(wrms-wmean*wmean)));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void addDataVariables(const_MatrixSliceRef data, VariableList &varList)
{
  try
  {
    varList.undefineVariable("index");
    for(UInt i=0; i<data.columns(); i++)
      addDataVariables("data"+i%"%i"s, data.column(i), varList);
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
    varList.setVariable("index", static_cast<Double>(row)); // index
    for(UInt i=0; i<data.columns(); i++)
      varList.setVariable("data"+i%"%i"s, data(row,i)); // "data(i)"
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
    varList.undefineVariable("index"); // index
    for(UInt i=0; i<data.columns(); i++)
      varList.undefineVariable("data"+i%"%i"s); // "data(i)"
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void addDataVariables(const GriddedData &grid, VariableList &varList)
{
  try
  {
    varList.undefineVariable("longitude");
    varList.undefineVariable("latitude");
    varList.undefineVariable("height");
    varList.undefineVariable("cartesianX");
    varList.undefineVariable("cartesianY");
    varList.undefineVariable("cartesianZ");
    varList.undefineVariable("area");
    varList.undefineVariable("index");
    for(UInt i=0; i<grid.values.size(); i++)
      addDataVariables("data"+i%"%i"s, Vector(grid.values.at(i)), Vector(grid.areas), varList);
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
    varList.setVariable("longitude", Double(L)*RAD2DEG);
    varList.setVariable("latitude", Double(B)*RAD2DEG);
    varList.setVariable("height", h);
    varList.setVariable("cartesianX", grid.points.at(row).x());
    varList.setVariable("cartesianY", grid.points.at(row).y());
    varList.setVariable("cartesianZ", grid.points.at(row).z());
    varList.setVariable("area", (grid.areas.size() ? grid.areas.at(row) : 0));
    varList.setVariable("index", static_cast<Double>(row));
    for(UInt i=0; i<grid.values.size(); i++)
      varList.setVariable("data"+i%"%i"s, grid.values.at(i).at(row));
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
    varList.undefineVariable("longitude");
    varList.undefineVariable("latitude");
    varList.undefineVariable("height");
    varList.undefineVariable("cartesianX");
    varList.undefineVariable("cartesianY");
    varList.undefineVariable("cartesianZ");
    varList.undefineVariable("area");
    varList.undefineVariable("index");
    for(UInt i=0; i<grid.values.size(); i++)
      varList.undefineVariable("data"+i%"%i"s);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void addDataVariables(const GriddedDataRectangular &grid, VariableList &varList)
{
  try
  {
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

    varList.undefineVariable("longitude");
    varList.undefineVariable("latitude");
    varList.undefineVariable("height");
    varList.undefineVariable("cartesianX");
    varList.undefineVariable("cartesianY");
    varList.undefineVariable("cartesianZ");
    varList.undefineVariable("index");
    for(UInt i=0; i<grid.values.size(); i++)
    {
      const std::string prefix = "data"+i%"%i"s;
      addDataVariables(prefix, grid.values.at(i), varList);

      Double wmean = 0;
      Double wrms  = 0;
      for(UInt z=0; z<rows; z++)
        for(UInt s=0; s<cols; s++)
        {
          const Double w = dLambda.at(s)*dPhi.at(z)/totalArea;
          wmean += w * grid.values.at(i)(z,s);
          wrms  += w * grid.values.at(i)(z,s) * grid.values.at(i)(z,s);
        }

      varList.setVariable(prefix+"wmean", wmean);
      varList.setVariable(prefix+"wrms",  std::sqrt(wrms));
      varList.setVariable(prefix+"wstd",  std::sqrt(rows*cols/(rows*cols-1.)*(wrms-wmean*wmean)));
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
    varList.setVariable("longitude", grid.longitudes.at(col)*RAD2DEG);
    varList.setVariable("latitude", grid.latitudes.at(row)*RAD2DEG);
    varList.setVariable("height", grid.heights.at(row));
    varList.setVariable("cartesianX", p.x());
    varList.setVariable("cartesianY", p.y());
    varList.setVariable("cartesianZ", p.z());
    varList.setVariable("index", static_cast<Double>(row*grid.longitudes.size()+col));
    for(UInt i=0; i<grid.values.size(); i++)
      varList.setVariable("data"+i%"%i"s, grid.values.at(i)(row,col)); // "data(i)"
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
    varList.undefineVariable("longitude");
    varList.undefineVariable("latitude");
    varList.undefineVariable("height");
    varList.undefineVariable("cartesianX");
    varList.undefineVariable("cartesianY");
    varList.undefineVariable("cartesianZ");
    varList.undefineVariable("index");
    for(UInt i=0; i<grid.values.size(); i++)
      varList.undefineVariable("data"+i%"%i"s); // "data(i)"
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/
