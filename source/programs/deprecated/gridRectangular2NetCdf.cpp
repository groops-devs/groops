/***********************************************/
/**
* @file gridRectangular2NetCdf.cpp
*
* @brief DEPRECATED. Please use GriddedData2NetCdf or GriddedDataTimeSeries2NetCdf instead.
*
* @author Andreas Kvas
* @author Torsten Mayer-Guerr
* @date 2018-06-17
*
* @deprecated Please use GriddedData2NetCdf or GriddedDataTimeSeries2NetCdf instead.
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
DEPRECATED. Please use \program{GriddedData2NetCdf} or \program{GriddedDataTimeSeries2NetCdf} instead.
)";

/***********************************************/

#include "programs/program.h"
#include "parser/dataVariables.h"
#include "inputOutput/fileNetCdf.h"
#include "files/fileGriddedData.h"
#include "classes/timeSeries/timeSeries.h"

/***** CLASS ***********************************/

/** @brief DEPRECATED. Please use GriddedData2NetCdf or GriddedDataTimeSeries2NetCdf instead.
* @ingroup programsConversionGroup */
class GridRectangular2NetCdf
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

GROOPS_REGISTER_PROGRAM(GridRectangular2NetCdf, SINGLEPROCESS, "DEPRECATED. Please use GriddedData2NetCdf or GriddedDataTimeSeries2NetCdf instead.", Deprecated)

/***********************************************/

template<> Bool readConfig(Config &config, const std::string &name, GridRectangular2NetCdf::Attribute &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  std::string choice;
  if(!readConfigChoice(config, name, choice, mustSet, defaultValue, annotation))
    return FALSE;
  if(readConfigChoiceElement(config, "text", choice, ""))
  {
    readConfig(config, "name",  var.name,  Config::MUSTSET, "", "");
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

template<> Bool readConfig(Config &config, const std::string &name, GridRectangular2NetCdf::DataVariable &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  if(!readConfigSequence(config, name, mustSet, defaultValue, annotation))
    return FALSE;
  readConfig(config, "selectDataField", var.dataField, Config::DEFAULT,  "0",     "input data column");
  readConfig(config, "name",            var.name,      Config::MUSTSET,  "data0", "netCDF variable name");
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

void GridRectangular2NetCdf::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName                  fileNameOut;
    std::vector<FileName>     fileNameIn;
    std::vector<DataVariable> dataVariables;
    std::vector<Attribute>    globalAttributes;
    TimeSeriesPtr             timeSeries;

    readConfig(config, "outputfileNetCdf",         fileNameOut,      Config::MUSTSET,  "", "file name of NetCDF output");
    readConfig(config, "inputfileGridRectangular", fileNameIn,       Config::MUSTSET,  "", "input grid sequence");
    readConfig(config, "times",                    timeSeries,       Config::DEFAULT,  "", "values for time axis (COARDS specification)");
    readConfig(config, "dataVariable",             dataVariables,    Config::MUSTSET,  "", "metadata for data variables");
    readConfig(config, "globalAttribute",          globalAttributes, Config::OPTIONAL, "", "additional meta data");
    if(isCreateSchema(config)) return;

    logWarning<<"DEPRECATED. Please use GriddedData2NetCdf or GriddedDataTimeSeries2NetCdf instead."<<Log::endl;

#ifdef GROOPS_DISABLE_NETCDF
    throw(Exception("Compiled without NetCDF library"));
#else
    // check input
    // -----------
    std::vector<Time> times = timeSeries->times();
    if(times.size() && (times.size() != fileNameIn.size()))
      throw(Exception("Mismatch in input file number and epoch count ("+fileNameIn.size()%"%i vs. "s+times.size()%"%i)."s));
    else if(!times.size() && (fileNameIn.size() > 1))
      throw(Exception("Specifying more than one input file without time axis results in undefined behavior."));

    // write attributes
    // ----------------
    NetCdf::OutFile file(fileNameOut);
    for(const auto &attr : globalAttributes)
    {
      if(attr.text.empty())
        file.addAttribute(attr.name, attr.dataType, attr.values);
      else
        file.addAttribute(attr.name, attr.text);
    }

    Bool firstGrid = TRUE;
    std::vector<Double> mjd;
    for(UInt idEpoch=0; idEpoch<fileNameIn.size(); idEpoch++)
    {
      GriddedDataRectangular grid;
      try
      {
        logStatus<<"read file <"<<fileNameIn.at(idEpoch)<<">"<<Log::endl;
        readFileGriddedData(fileNameIn.at(idEpoch), grid);
      }
      catch(std::exception &e)
      {
        logWarning<<e.what()<<Log::endl;
        continue;
      }

      // create dimensions and variables
      // -------------------------------
      if(firstGrid)
      {
        firstGrid = FALSE;
        if(times.size())
        {
          auto dim = file.addDimension("time");
          auto var = file.addVariable("time", NetCdf::DOUBLE, {dim});
          var.addAttribute("units", "days since 1858-11-17 00:00:00"); // MJD
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
        }
      } // if(firstGrid)

      for(auto &data : dataVariables)
      {
        Vector values = flatten(grid.values.at(data.dataField).trans());
        if(times.size())
          data.variable.setValues({mjd.size(), 0, 0}, {1, grid.latitudes.size(), grid.longitudes.size()}, values);
        else
          data.variable.setValues({0, 0}, {grid.latitudes.size(), grid.longitudes.size()}, values);
      }
      if(times.size())
        mjd.push_back(times.at(idEpoch).mjd());
    }

    if(mjd.size())
      file.variable("time").setValues(mjd);
#endif
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

