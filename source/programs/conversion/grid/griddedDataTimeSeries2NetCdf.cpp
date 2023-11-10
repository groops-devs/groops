/***********************************************/
/**
* @file griddedDataTimeSeries2NetCdf.cpp
*
* @brief Converts a griddedDataTimeSeries file to a netCDF file.
*
* @author Andreas Kvas
* @author Torsten Mayer-Guerr
* @date 2018-06-17
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Read a \configFile{inputfileGriddedDataTimeSeries}{griddedDataTimeSeries}
and converts it to a COARDS compliant NetCDF file.

The output data can be defined with \config{dataVariable}.
You should add at least the attributes \verb|units|, \verb|long_name|, and maybe \verb|_FillValue|
to the variables. The \config{dataVariable:inputColumn} selects the data from the input file.

If \configClass{timeSeries}{timeSeriesType} is not set
the temporal nodal points from the inputfile are used.

See also \program{NetCdfInfo}, \program{GriddedData2NetCdf}, \program{NetCdf2GriddedDataTimeSeries}.
)";

/***********************************************/

#include "programs/program.h"
#include "inputOutput/fileNetCdf.h"
#include "files/fileGriddedDataTimeSeries.h"
#include "classes/timeSeries/timeSeries.h"
#include "misc/miscGriddedData.h"

/***** CLASS ***********************************/

/** @brief Converts a griddedDataTimeSeries file to a netCDF file.
* @ingroup programsConversionGroup */
class GriddedDataTimeSeries2NetCdf
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
    UInt                   dataField;
    NetCdf::DataType       dataType;
    std::vector<Attribute> attributes;
    NetCdf::Variable       variable;
  };

  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GriddedDataTimeSeries2NetCdf, SINGLEPROCESS, "Converts a griddedDataTimeSeries file to a netCDF file", Conversion, Grid)

/***********************************************/

template<> Bool readConfig(Config &config, const std::string &name, GriddedDataTimeSeries2NetCdf::Attribute &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
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

template<> Bool readConfig(Config &config, const std::string &name, GriddedDataTimeSeries2NetCdf::DataVariable &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  if(!readConfigSequence(config, name, mustSet, defaultValue, annotation))
    return FALSE;
  readConfig(config, "name",        var.name,      Config::MUSTSET,  "data0", "netCDF variable name");
  readConfig(config, "inputColumn", var.dataField, Config::DEFAULT,  "0",     "input data column");
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

void GriddedDataTimeSeries2NetCdf::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName                  fileNameOut, fileNameIn;
    std::vector<DataVariable> dataVariables;
    std::vector<Attribute>    globalAttributes;
    TimeSeriesPtr             timeSeries;

    readConfig(config, "outputfileNetCdf",               fileNameOut,      Config::MUSTSET,  "", "file name of NetCDF output");
    readConfig(config, "inputfileGriddedDataTimeSeries", fileNameIn,       Config::MUSTSET,  "", "");
    readConfig(config, "timeSeries",                     timeSeries,       Config::OPTIONAL, "", "otherwise times from inputfile are used");
    readConfig(config, "dataVariable",                   dataVariables,    Config::MUSTSET,  "", "metadata for data variables");
    readConfig(config, "globalAttribute",                globalAttributes, Config::OPTIONAL, "", "additional meta data");
    if(isCreateSchema(config)) return;

#ifdef GROOPS_DISABLE_NETCDF
    throw(Exception("Compiled without NetCDF library"));
#else

    logStatus<<"read gridded data time series <"<<fileNameIn<<">"<<Log::endl;
    InFileGriddedDataTimeSeries fileGriddedDataTimeSeries(fileNameIn);
    GriddedDataRectangular grid;
    if(!grid.init(fileGriddedDataTimeSeries.grid()))
      throw(Exception("GriddedData must be a rectangle grid"));
    MiscGriddedData::printStatistics(grid);
    std::vector<Time> times = fileGriddedDataTimeSeries.times();
    if(timeSeries)
      times = timeSeries->times();

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
      NetCdf::Dimension dim = file.addDimension("time"); // NC_UNLIMITED
      NetCdf::Variable  var = file.addVariable("time", NetCdf::DOUBLE, {dim});
      var.addAttribute("units", "days since 1858-11-17 00:00:00"); // MJD
    }
    {
      NetCdf::Dimension dim = file.addDimension("lat", grid.latitudes.size());
      NetCdf::Variable  var = file.addVariable("lat", NetCdf::DOUBLE, {dim});
      var.addAttribute("units", "degrees_north");
      Vector lats(grid.latitudes.size());
      for(UInt i=0; i<grid.latitudes.size(); i++)
        lats(i) = RAD2DEG * grid.latitudes.at(i);
      var.setValues(lats);
    }
    {
      NetCdf::Dimension dim = file.addDimension("lon", grid.longitudes.size());
      NetCdf::Variable  var = file.addVariable("lon", NetCdf::DOUBLE, {dim});
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
    }

    Single::forEach(times.size(), [&](UInt idEpoch)
    {
      Matrix values = fileGriddedDataTimeSeries.data(times.at(idEpoch));
      for(auto &data : dataVariables)
        data.variable.setValues({idEpoch, 0, 0}, {1, grid.latitudes.size(), grid.longitudes.size()}, Vector(values.column(data.dataField)));
    });

    Vector mjd(times.size());
    for(UInt i=0; i<times.size(); i++)
      mjd(i) = times.at(i).mjd();
    file.variable("time").setValues(mjd);
#endif
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

