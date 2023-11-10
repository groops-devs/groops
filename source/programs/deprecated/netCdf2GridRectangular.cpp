/***********************************************/
/**
* @file netCdf2GridRectangular.cpp
*
* @brief DEPRECATED. Please use NetCdf2GriddedData or NetCdf2GriddedDataTimeSeries instead.
*
* @author Andreas Kvas
* @author Torsten Mayer-Guerr
* @date 2018-06-17
*
* @deprecated Please use NetCdf2GriddedData or NetCdf2GriddedDataTimeSeries instead.
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
DEPRECATED. Please use \program{NetCdf2GriddedData} or \program{NetCdf2GriddedDataTimeSeries} instead.
)";

/***********************************************/

#include "programs/program.h"
#include "inputOutput/fileNetCdf.h"
#include "files/fileGriddedData.h"

/***** CLASS ***********************************/

/** @brief DEPRECATED. Please use NetCdf2GriddedData or NetCdf2GriddedDataTimeSeries instead.
* @ingroup programsConversionGroup */
class NetCdf2GridRectangular
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(NetCdf2GridRectangular, SINGLEPROCESS, "DEPRECATED. Please use NetCdf2GriddedData or NetCdf2GriddedDataTimeSeries instead.", Deprecated)

/***********************************************/

void NetCdf2GridRectangular::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName    outName, inName;
    std::string loopVar;
    std::string lonName, latName, timeName;
    std::vector<std::string> dataName;
    Double      a, f;

    readConfig(config, "outputfileGridRectangular", outName,  Config::MUSTSET,  "",         "One grid for each epoch in the NetCDF file is written. Use loopTimeVariable as template.");
    readConfig(config, "loopTimeVariable",          loopVar,  Config::MUSTSET,  "loopTime", "");
    readConfig(config, "inputfileNetCdf",           inName,   Config::MUSTSET,  "",         "");
    readConfig(config, "variableNameLongitude",     lonName,  Config::MUSTSET,  "lon",      "name of NetCDF variable");
    readConfig(config, "variableNameLatitude",      latName,  Config::MUSTSET,  "lat",      "name of NetCDF variable");
    readConfig(config, "variableNameTime",          timeName, Config::OPTIONAL, "time",     "name of NetCDF variable (leave blank for static grids)");
    readConfig(config, "variableNameData",          dataName, Config::OPTIONAL, "",         "name of NetCDF variable");
    readConfig(config, "R",                         a,        Config::DEFAULT,  STRING_DEFAULT_GRS80_a, "reference radius for ellipsoidal coordinates");
    readConfig(config, "inverseFlattening",         f,        Config::DEFAULT,  STRING_DEFAULT_GRS80_f, "reference flattening for ellipsoidal coordinates");
    if(isCreateSchema(config)) return;

    logWarning<<"DEPRECATED. Please use NetCdf2GriddedData or NetCdf2GriddedDataTimeSeries instead."<<Log::endl;

#ifdef GROOPS_DISABLE_NETCDF
    throw(Exception("Compiled without NetCDF library"));
#else
    // add variable for output files
    // -----------------------------
    VariableList fileNameVariableList;
    fileNameVariableList.undefineVariable(loopVar);

    // open netCDF file
    // ----------------
    logStatus<<"read file <"<<inName<<">"<<Log::endl;
    NetCdf::InFile file(inName);
    NetCdf::Variable  lon    = file.variable(lonName);
    NetCdf::Variable  lat    = file.variable(latName);
    NetCdf::Dimension dimLon = lon.dimensions().at(0);
    NetCdf::Dimension dimLat = lat.dimensions().at(0);

    // set up grid
    // -----------
    GriddedDataRectangular grid;
    grid.ellipsoid  = Ellipsoid(a, f);
    grid.longitudes = NetCdf::convertAngles(lon.values());
    grid.latitudes  = NetCdf::convertAngles(lat.values());
    grid.heights.resize(grid.latitudes.size(), 0.0);

    // set up time axis
    // ----------------
    std::vector<Time> epochs(1);
    NetCdf::Dimension dimTime;
    if(!timeName.empty())
    {
      auto var = file.variable(timeName);
      dimTime  = var.dimensions().at(0);
      epochs   = NetCdf::convertTimes(var.values(), var.attribute("units").value());
    }

    // data variables
    // --------------
    Single::forEach(epochs.size(), [&](UInt idEpoch)
    {
      fileNameVariableList.setVariable(loopVar, epochs.at(idEpoch).mjd());

      grid.values.clear();
      for(const std::string &name : dataName)
      {
        auto var  = file.variable(name);
        auto dims = var.dimensions();

        std::vector<UInt> start(dims.size(), 0);
        std::vector<UInt> count;
        for(auto &dim : dims)
          count.push_back(dim.length());

        if((dims.size() != 2) && (timeName.empty() || (dims.size() != 3)))
           throw(Exception("variable <"+name+"> has wrong dimensions"));

        if(!timeName.empty() && (dims.size() > 2))
        {
          if(dims.at(0) != dimTime)
            throw(Exception("variable <"+name+"> must have time as first dimension"));
          start.at(0) = idEpoch;
          count.at(0) = 1;
        }

        Matrix values = reshape(var.values(start, count), count.at(dims.size()-1), count.at(dims.size()-2));
        if((dims.at(dims.size()-1) == dimLat) && (dims.at(dims.size()-2) == dimLon))
          grid.values.push_back(values);
        else if((dims.at(dims.size()-1) == dimLon) && (dims.at(dims.size()-2) == dimLat))
          grid.values.push_back(values.trans());
        else
           throw(Exception("variable <"+name+"> must have ("+latName+", "+lonName+") dimensions"));
      }

      writeFileGriddedData(outName(fileNameVariableList), grid);
    });
#endif
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
