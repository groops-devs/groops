/***********************************************/
/**
* @file griddedData2NetCdf.cpp
*
* @brief Convert gridded data to a netCDF file.
*
* @author Andreas Kvas
* @author Torsten Mayer-Guerr
* @date 2018-06-17
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program converts a \configFile{inputfileGriddedData}{griddedData}
to a COARDS compliant NetCDF file. The output data can be defined with \config{dataVariable}.
You should add at least the attributes \verb|units|, \verb|long_name|, and maybe \verb|_FillValue|
to the variables. For the \config{dataVariable:value} the standard \reference{dataVariables}{general.parser:dataVariables}
are available to select the data columns of \configFile{inputfileGriddedData}{griddedData}.

See also \program{NetCdfInfo}, \program{GriddedDataTimeSeries2NetCdf}, \program{NetCdf2GriddedData}.
)";

/***********************************************/

#include "programs/program.h"
#include "parser/dataVariables.h"
#include "inputOutput/fileNetCdf.h"
#include "files/fileGriddedData.h"

/***** CLASS ***********************************/

/** @brief Convert gridded data to a netCDF file.
* @ingroup programsConversionGroup */
class GriddedData2NetCdf
{
public:
  class Attribute
  {
  public:
    std::string         name;
    std::string         text;
    std::vector<Double> values;
    NetCdf::DataType    dataType;
  };

  class DataVariable
  {
  public:
    std::string            name;
    ExpressionVariablePtr  valueExpr;
    Vector                 values;
    NetCdf::DataType       dataType;
    std::vector<Attribute> attributes;
    NetCdf::Variable       variable;
  };

  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GriddedData2NetCdf, SINGLEPROCESS, "Convert gridded data to a netCDF file", Conversion, Grid)

/***********************************************/

template<> Bool readConfig(Config &config, const std::string &name, GriddedData2NetCdf::Attribute &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  std::string choice;
  if(!readConfigChoice(config, name, choice, mustSet, defaultValue, annotation))
    return FALSE;
  if(readConfigChoiceElement(config, "text", choice, ""))
  {
    readConfig(config, "name",  var.name, Config::MUSTSET, "", "");
    readConfig(config, "value", var.text, Config::MUSTSET, "", "");
  }
  if(readConfigChoiceElement(config, "value",  choice, ""))
  {
    std::string choice;
    readConfig(config, "name",  var.name,   Config::MUSTSET, "", "");
    readConfig(config, "value", var.values, Config::MUSTSET, "", "");
    if(readConfigChoice(config, "dataType", choice, Config::MUSTSET, "", ""))
    {
      if(readConfigChoiceElement(config, "double", choice, "")) var.dataType = NetCdf::DOUBLE;
      if(readConfigChoiceElement(config, "float",  choice, "")) var.dataType = NetCdf::FLOAT;
      if(readConfigChoiceElement(config, "int",    choice, "")) var.dataType = NetCdf::INT;
      endChoice(config);
    }
  }
  endChoice(config);
  return TRUE;
}

/***********************************************/

template<> Bool readConfig(Config &config, const std::string &name, GriddedData2NetCdf::DataVariable &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  if(!readConfigSequence(config, name, mustSet, defaultValue, annotation))
    return FALSE;
  readConfig(config, "name",  var.name,      Config::MUSTSET,  "data0", "netCDF variable name");
  readConfig(config, "value", var.valueExpr, Config::DEFAULT,  "data0", "expression (variables 'height', 'data', 'longitude', 'latitude' and, 'area' are taken from the gridded data");
  std::string choice;
  if(readConfigChoice(config, "dataType", choice, Config::MUSTSET, "", ""))
  {
    if(readConfigChoiceElement(config, "double", choice, "")) var.dataType = NetCdf::DOUBLE;
    if(readConfigChoiceElement(config, "float",  choice, "")) var.dataType = NetCdf::FLOAT;
    if(readConfigChoiceElement(config, "int",    choice, "")) var.dataType = NetCdf::INT;
    endChoice(config);
  }
  readConfig(config, "attribute", var.attributes, Config::OPTIONAL, "", "netCDF attributes");
  endSequence(config);
  return TRUE;
}

/***********************************************/

void GriddedData2NetCdf::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName                  fileNameOut, fileNameIn;
    std::vector<DataVariable> dataVariables;
    std::vector<Attribute>    globalAttributes;

    readConfig(config, "outputfileNetCdf",     fileNameOut,      Config::MUSTSET,  "", "file name of NetCDF output");
    readConfig(config, "inputfileGriddedData", fileNameIn,       Config::MUSTSET,  "", "input grid sequence");
    readConfig(config, "dataVariable",         dataVariables,    Config::MUSTSET,  "", "metadata for data variables");
    readConfig(config, "globalAttribute",      globalAttributes, Config::OPTIONAL, "", "additional meta data");
    if(isCreateSchema(config)) return;

#ifdef GROOPS_DISABLE_NETCDF
    throw(Exception("Compiled without NetCDF library"));
#else

    // read grid
    // ---------
    logStatus<<"read gridded data <"<<fileNameIn<<">"<<Log::endl;
    GriddedDataRectangular grid;
    readFileGriddedData(fileNameIn, grid);

    // create data variables
    // ---------------------
    VariableList varList;
    addDataVariables(grid, varList);
    for(auto &data : dataVariables)
      data.valueExpr->simplify(varList);

    for(auto &data : dataVariables)
      data.values = Vector(grid.longitudes.size() * grid.latitudes.size());

    UInt idx = 0;
    for(UInt i=0; i<grid.latitudes.size(); i++)
      for(UInt k=0; k<grid.longitudes.size(); k++)
      {
        evaluateDataVariables(grid, i, k, varList);
        for(auto &data : dataVariables)
          data.values(idx) = data.valueExpr->evaluate(varList);
        idx++;
      }

    // write attributes
    // ----------------
    logStatus<<"write NETCDF file <"<<fileNameOut<<">"<<Log::endl;
    NetCdf::OutFile file(fileNameOut);
    for(const auto &attr : globalAttributes)
    {
      if(attr.text.empty())
        file.addAttribute(attr.name, attr.dataType, attr.values);
      else
        file.addAttribute(attr.name, attr.text);
    }
    {
      auto dim = file.addDimension("lat", grid.latitudes.size());
      auto var = file.addVariable("lat", NetCdf::DOUBLE, {dim});
      var.addAttribute("units", "degrees_north");
      Vector lats(grid.latitudes.size());
      for(UInt i=0; i<grid.latitudes.size(); i++)
        lats(i) = RAD2DEG * grid.latitudes.at(i);
      var.setValues(lats);
    }
    {
      auto dim = file.addDimension("lon", grid.longitudes.size());
      auto var = file.addVariable("lon", NetCdf::DOUBLE, {dim});
      var.addAttribute("units", "degrees_east");
      Vector lons(grid.longitudes.size());
      for(UInt i=0; i<grid.longitudes.size(); i++)
        lons(i) = RAD2DEG * grid.longitudes.at(i);
      var.setValues(lons);
    }
    for(auto &data : dataVariables)
    {
      data.variable = file.addVariable(data.name, data.dataType, file.dimensions());
      for(auto &attr : data.attributes)
      {
        if(attr.text.empty())
          data.variable.addAttribute(attr.name, attr.dataType, attr.values);
        else
          data.variable.addAttribute(attr.name, attr.text);
      }
      data.variable.setValues({0, 0}, {grid.latitudes.size(), grid.longitudes.size()}, data.values);
    }
#endif
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

