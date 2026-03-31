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
        if((std::find(dims.begin(), dims.end(), dimLat) == dims.end()) ||
           (std::find(dims.begin(), dims.end(), dimLon) == dims.end()))
          continue;
        UInt countDims = 0;
        for(auto &dim : dims)
          if((dim != dimLat) && (dim != dimLon) && (!timeName.empty() && (dim != dimTime)) && (dim.length() > 1))
            countDims++;
        if(countDims <= 1)
          dataNames.push_back(var.name());
      }
      if(!dataNames.size())
        logWarning<<"No suitable variables found"<<Log::endl;
    }

    // read data variables
    // -------------------
    UInt idxCol = 0;
    for(const std::string &name : dataNames)
    {
      NetCdf::Variable var = file.variable(name);
      std::vector<NetCdf::Dimension> dims = var.dimensions();

      std::stringstream ss;
      if(dims.size())
      {
        ss<<"("<<dims.at(0).name();
        for(UInt i=1; i<dims.size(); i++)
          ss<<", "<<dims.at(i).name();
        ss<<")";
      }
      logInfo<<"  data"<<idxCol<<" = "<<var.name()<<ss.str()<<Log::endl;
      for(auto &attr : var.attributes())
        logInfo<<"    - "<<attr.name()<<" value = "<<attr.value()<<Log::endl;

      if((std::find(dims.begin(), dims.end(), dimLat) == dims.end()) ||
         (std::find(dims.begin(), dims.end(), dimLon) == dims.end()))
        throw(Exception("variable <"+var.name()+"> must have at least lat, lon dimensions"));

      std::vector<UInt> start(dims.size(), 0);
      std::vector<UInt> count(dims.size(), 0);
      UInt incLat = 0, incLon = 0, incCol = 0;
      UInt countColumns = 1;
      UInt inc = 1;
      for(UInt k=dims.size(); k-->0;)
      {
        count.at(k) = dims.at(k).length();
        if(!timeName.empty() && (dims.at(k) == dimTime))
        {
          start.at(k) = idEpoch;
          count.at(k) = 1;
        }
        else
        {
          if(dims.at(k) == dimLat)
            incLat = inc;
          else if(dims.at(k) == dimLon)
            incLon = inc;
          else if(dims.at(k).length() > 1)
          {
            if(countColumns > 1)
              throw(Exception("variable <"+var.name()+"> has wrong dimensions"));
            incCol = inc;
            countColumns = dims.at(k).length();
          }
          inc *= dims.at(k).length();
        }
      }
      idxCol += countColumns;

      // get data
      const Vector values = var.values(start, count);
      for(UInt j=0; j<countColumns; j++)
      {
        Matrix data(grid.latitudes.size(), grid.longitudes.size());
        for(UInt i=0; i<grid.latitudes.size(); i++)
          for(UInt k=0; k<grid.longitudes.size(); k++)
            data(i, k) = values(incLat*i+incLon*k+incCol*j);
        grid.values.push_back(data);
      }
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
