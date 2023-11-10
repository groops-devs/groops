/***********************************************/
/**
* @file netCdf2GriddedData.cpp
*
* @brief Convert COARDS compliant grids to gridded data.
*
* @author Andreas Kvas
* @author Torsten Mayer-Guerr
* @date 2023-07-04
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program converts a COARDS compliant NetCDF file into an
\configFile{outputfileGriddedData}{griddedData}.
If no specific input \config{variableNameData} are selected all suitable data are used.

If the NETCDF file contains a time axis (\config{variableNameData}) an specific epoch
can be selected with \config{time}. The nearest epoch in file is used.

See also \program{NetCdfInfo}, \program{GriddedData2NetCdf}, \program{NetCdf2GriddedDataTimeSeries}.
)";

/***********************************************/

#include "programs/program.h"
#include "inputOutput/fileNetCdf.h"
#include "files/fileGriddedData.h"
#include "misc/miscGriddedData.h"

/***** CLASS ***********************************/

/** @brief Convert COARDS compliant grids to gridded data.
* @ingroup programsConversionGroup */
class NetCdf2GriddedData
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(NetCdf2GriddedData, SINGLEPROCESS, "Convert COARDS compliant grids to gridded data", Conversion, Grid)

/***********************************************/

void NetCdf2GriddedData::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName    fileNameGriddedData, fileNameNetcdf;
    std::string lonName, latName, timeName;
    std::vector<std::string> dataNames;
    Time        time;
    Double      a, f;

    readConfig(config, "outputfileGriddedData", fileNameGriddedData, Config::MUSTSET,  "",     "");
    readConfig(config, "inputfileNetCdf",       fileNameNetcdf,      Config::MUSTSET,  "",     "");
    readConfig(config, "variableNameLongitude", lonName,             Config::MUSTSET,  "lon",  "name of NetCDF variable");
    readConfig(config, "variableNameLatitude",  latName,             Config::MUSTSET,  "lat",  "name of NetCDF variable");
    readConfig(config, "variableNameTime",      timeName,            Config::OPTIONAL, "time", "if with time axis: name of NetCDF variable");
    readConfig(config, "variableNameData",      dataNames,           Config::OPTIONAL, "",     "data variables, otherwise all suitable data are used");
    readConfig(config, "time",                  time,                Config::OPTIONAL, "",     "if with time axis: nearest epoch is used");
    readConfig(config, "R",                     a,                   Config::DEFAULT,  STRING_DEFAULT_GRS80_a, "reference radius for ellipsoidal coordinates");
    readConfig(config, "inverseFlattening",     f,                   Config::DEFAULT,  STRING_DEFAULT_GRS80_f, "reference flattening for ellipsoidal coordinates");
    if(isCreateSchema(config)) return;

#ifdef GROOPS_DISABLE_NETCDF
    throw(Exception("Compiled without NetCDF library"));
#else
    // open netCDF file
    // ----------------
    logStatus<<"read NETCDF file <"<<fileNameNetcdf<<">"<<Log::endl;
    NetCdf::InFile file(fileNameNetcdf);
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
    UInt dimSize = 2; // lat, lon
    UInt idEpoch = 0;
    NetCdf::Dimension dimTime;
    if(!timeName.empty() && file.hasVariable(timeName))
    {
      auto var = file.variable(timeName);
      dimTime  = var.dimensions().at(0);
      std::vector<Time> epochs = NetCdf::convertTimes(var.values(), var.attribute("units").value());
      idEpoch = std::distance(epochs.begin(), std::min_element(epochs.begin(), epochs.end(),
                              [&](Time &t1, Time &t2) {return std::fabs((t1-time).mjd()) < std::fabs((t2-time).mjd());}));
      logInfo<<"  epoch used: "<<epochs.at(idEpoch).dateTimeStr()<<Log::endl;
      dimSize++; // time, lat, lon
    }
    else
      timeName = std::string();

    // set up data columns
    // -------------------
    if(!dataNames.size())
    {
      for(auto &var : file.variables())
      {
        auto dims = var.dimensions();
        if((dims.size() != dimSize) ||
           (std::find(dims.begin(), dims.end(), dimLat) == dims.end()) ||
           (std::find(dims.begin(), dims.end(), dimLon) == dims.end()) ||
           (!timeName.empty() && (std::find(dims.begin(), dims.end(), dimTime) == dims.end())))
          continue;
        dataNames.push_back(var.name());
      }
      if(!dataNames.size())
        logWarning<<"No suitable variables found"<<Log::endl;
    }


    // read data variables
    // -------------------
    UInt idx = 0;
    for(const std::string &name : dataNames)
    {
      auto var  = file.variable(name);
      auto dims = var.dimensions();

      logInfo<<"  data"<<idx++<<" = "<<var.name()<<Log::endl;
      for(auto &attr : var.attributes())
        logInfo<<"    - "<<attr.name()<<" value = "<<attr.value()<<Log::endl;

      std::vector<UInt> start(dims.size(), 0);
      std::vector<UInt> count;
      for(auto &dim : dims)
        count.push_back(dim.length());

      if(dims.size() != dimSize)
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

    // write grid
    // ----------
    logStatus<<"write gridded data to <"<<fileNameGriddedData<<">"<<Log::endl;
    writeFileGriddedData(fileNameGriddedData, grid);
    MiscGriddedData::printStatistics(grid);
#endif
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
