/***********************************************/
/**
* @file netCdf2GriddedDataTimeSeries.cpp
*
* @brief Convert a NetCDF file to a GriddedDataTimeSeries file.
*
* @author Andreas Kvas
* @author Torsten Mayer-Guerr
* @date 2023-07-04
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program converts a COARDS compliant NetCDF file into
\configFile{outputfileGriddedDataTimeSeries}{griddedDataTimeSeries}.
If no specific input \config{variableNameData} are selected all suitable data are used.

See also \program{NetCdfInfo}, \program{NetCdf2GriddedData}, \program{GriddedDataTimeSeries2NetCdf}.
)";

/***********************************************/

#include "programs/program.h"
#include "inputOutput/fileNetCdf.h"
#include "files/fileGriddedDataTimeSeries.h"
#include "misc/miscGriddedData.h"

/***** CLASS ***********************************/

/** @brief Convert a NetCDF file to a GriddedDataTimeSeries file.
* @ingroup programsConversionGroup */
class NetCdf2GriddedDataTimeSeries
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(NetCdf2GriddedDataTimeSeries, SINGLEPROCESS, "Convert a NetCDF file to a GriddedDataTimeSeries file", Conversion, Grid)

/***********************************************/

void NetCdf2GriddedDataTimeSeries::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName    fileNameOut, fileNameNetcdf;
    std::string lonName, latName, timeName;
    std::vector<std::string> dataNames;
    Double      a, f;

    readConfig(config, "outputfileGriddedDataTimeSeries", fileNameOut,    Config::MUSTSET,  "",     "");
    readConfig(config, "inputfileNetCdf",                 fileNameNetcdf, Config::MUSTSET,  "",     "");
    readConfig(config, "variableNameLongitude",           lonName,        Config::MUSTSET,  "lon",  "name of NetCDF variable");
    readConfig(config, "variableNameLatitude",            latName,        Config::MUSTSET,  "lat",  "name of NetCDF variable");
    readConfig(config, "variableNameTime",                timeName,       Config::MUSTSET,  "time", "name of NetCDF variable)");
    readConfig(config, "variableNameData",                dataNames,      Config::OPTIONAL, "",     "data variables, otherwise all suitable data are used");
    readConfig(config, "R",                               a,              Config::DEFAULT,  STRING_DEFAULT_GRS80_a, "reference radius for ellipsoidal coordinates");
    readConfig(config, "inverseFlattening",               f,              Config::DEFAULT,  STRING_DEFAULT_GRS80_f, "reference flattening for ellipsoidal coordinates");
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
    GriddedData grid;
    GriddedDataRectangular griddedDataRectangular;
    griddedDataRectangular.ellipsoid  = Ellipsoid(a, f);
    griddedDataRectangular.longitudes = NetCdf::convertAngles(lon.values());
    griddedDataRectangular.latitudes  = NetCdf::convertAngles(lat.values());
    griddedDataRectangular.heights.resize(griddedDataRectangular.latitudes.size(), 0.0);
    griddedDataRectangular.convert(grid);
    MiscGriddedData::printStatistics(grid);

    // set up time axis
    // ----------------
    NetCdf::Variable  time    = file.variable(timeName);
    NetCdf::Dimension dimTime = time.dimensions().at(0);
    std::vector<Time> times   = NetCdf::convertTimes(time.values(), time.attribute("units").value());

    // set up data columns
    // -------------------
    if(!dataNames.size())
    {
      for(auto &var : file.variables())
      {
        auto dims = var.dimensions();
        if((dims.size() != 3) ||
           (std::find(dims.begin(), dims.end(), dimLat) == dims.end()) ||
           (std::find(dims.begin(), dims.end(), dimLon) == dims.end()) ||
           (std::find(dims.begin(), dims.end(), dimTime) == dims.end()))
          continue;
        dataNames.push_back(var.name());
      }
      if(!dataNames.size())
        logWarning<<"No suitable variables found"<<Log::endl;
    }

    std::vector<NetCdf::Variable>  vars(dataNames.size());
    std::vector<std::vector<NetCdf::Dimension>> dims(dataNames.size());
    UInt idx = 0;
    for(UInt i=0; i<vars.size(); i++)
    {
      vars.at(i) = file.variable(dataNames.at(i));
      dims.at(i) = vars.at(i).dimensions();
      logInfo<<"  data"<<idx++<<" = "<<vars.at(i).name()<<Log::endl;
      for(auto &attr : vars.at(i).attributes())
        logInfo<<"    - "<<attr.name()<<" value = "<<attr.value()<<Log::endl;
      if(dims.at(i).size() != 3)
        throw(Exception("variable <"+vars.at(i).name()+"> has wrong dimensions"));
      if(dims.at(i).at(0) != dimTime)
        throw(Exception("variable <"+vars.at(i).name()+"> must have time as first dimension"));
    }

    // read data variables
    // -------------------
    std::vector<Matrix> data(times.size(), Matrix(griddedDataRectangular.longitudes.size()*griddedDataRectangular.latitudes.size(), dataNames.size()));
    Single::forEach(times.size(), [&](UInt idEpoch)
    {
      for(UInt i=0; i<vars.size(); i++)
      {
        std::vector<UInt> start(dims.at(i).size(), 0);
        std::vector<UInt> count;
        for(auto &dim : dims.at(i))
          count.push_back(dim.length());
        start.at(0) = idEpoch;
        count.at(0) = 1;

        if((dims.at(i).at(1) == dimLon) && (dims.at(i).at(2) == dimLat))
          copy(vars.at(i).values(start, count), data.at(idEpoch).column(i));
        else if((dims.at(i).at(1) == dimLat) && (dims.at(i).at(2) == dimLon))
          reshape(reshape(vars.at(i).values(start, count), count.at(2), count.at(1)), data.at(idEpoch).column(i));
        else
          throw(Exception("variable <"+vars.at(i).name()+"> must have ("+latName+", "+lonName+") dimensions"));
      }
    });

    // write grid
    // ----------
    logStatus<<"write gridded data time series to <"<<fileNameOut<<">"<<Log::endl;
    writeFileGriddedDataTimeSeries(fileNameOut, 1/*splineDegree*/, times, grid, data);
#endif
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
