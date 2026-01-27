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
    GriddedDataRectangular griddedDataRectangular;
    griddedDataRectangular.ellipsoid  = Ellipsoid(a, f);
    griddedDataRectangular.longitudes = NetCdf::convertAngles(lon.values());
    griddedDataRectangular.latitudes  = NetCdf::convertAngles(lat.values());
    griddedDataRectangular.heights.resize(griddedDataRectangular.latitudes.size(), 0.0);
    MiscGriddedData::printStatistics(griddedDataRectangular);
    GriddedData grid(griddedDataRectangular);

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
        if((std::find(dims.begin(), dims.end(), dimLat)  == dims.end()) ||
           (std::find(dims.begin(), dims.end(), dimLon)  == dims.end()) ||
           (std::find(dims.begin(), dims.end(), dimTime) == dims.end()))
          continue;
        UInt countDims = 0;
        for(auto &dim : dims)
          if((dim != dimLat) && (dim != dimLon) && (dim != dimTime) && (dim.length() > 1))
            countDims++;
        if(countDims <= 1)
          dataNames.push_back(var.name());
      }
      if(!dataNames.size())
        logWarning<<"No suitable variables found"<<Log::endl;
    }

    class Var
    {
    public:
      NetCdf::Variable  var;
      std::vector<UInt> start, count;
      UInt              idxTime;
      UInt              incLat, incLon, incCol;
      UInt              countColumns;
    };

    std::vector<Var> vars(dataNames.size());
    UInt countColumns = 0;
    for(UInt i=0; i<vars.size(); i++)
    {
      vars.at(i).var = file.variable(dataNames.at(i));
      std::vector<NetCdf::Dimension> dims = vars.at(i).var.dimensions();

      std::stringstream ss;
      if(dims.size())
      {
        ss<<"("<<dims.at(0).name();
        for(UInt i=1; i<dims.size(); i++)
          ss<<", "<<dims.at(i).name();
        ss<<")";
      }
      logInfo<<"  data"<<countColumns<<" = "<<vars.at(i).var.name()<<ss.str()<<Log::endl;
      for(auto &attr : vars.at(i).var.attributes())
        logInfo<<"    - "<<attr.name()<<" value = "<<attr.value()<<Log::endl;

      if((std::find(dims.begin(), dims.end(), dimLat)  == dims.end()) ||
         (std::find(dims.begin(), dims.end(), dimLon)  == dims.end()) ||
         (std::find(dims.begin(), dims.end(), dimTime) == dims.end()))
        throw(Exception("variable <"+vars.at(i).var.name()+"> must have at least time, lat, lon dimensions"));

      vars.at(i).start = vars.at(i).count = std::vector<UInt>(dims.size(), 0);
      vars.at(i).countColumns = 1;
      UInt inc = 1;
      for(UInt k=dims.size(); k-->0;)
      {
        vars.at(i).count.at(k) = dims.at(k).length();
        if(dims.at(k) == dimTime)
        {
          vars.at(i).idxTime = k;
        }
        else
        {
          if(dims.at(k) == dimLat)
            vars.at(i).incLat = inc;
          else if(dims.at(k) == dimLon)
            vars.at(i).incLon = inc;
          else if(dims.at(k).length() > 1)
          {
            if(countColumns > 1)
              throw(Exception("variable <"+vars.at(i).var.name()+"> has wrong dimensions"));
            vars.at(i).incCol = inc;
            vars.at(i).countColumns = dims.at(k).length();
          }
          inc *= dims.at(k).length();
        }
      }
      countColumns += vars.at(i).countColumns;
    }

    // read data variables
    // -------------------
    std::vector<Matrix> data(times.size(), Matrix(griddedDataRectangular.longitudes.size()*griddedDataRectangular.latitudes.size(), countColumns));
    Single::forEach(times.size(), [&](UInt idEpoch)
    {
      UInt idxCol = 0;
      for(auto &var : vars)
      {
        var.start.at(var.idxTime) = idEpoch;
        var.count.at(var.idxTime) = 1;
        const Vector values = var.var.values(var.start, var.count);
        for(UInt i=0; i<griddedDataRectangular.latitudes.size(); i++)
          for(UInt k=0; k<griddedDataRectangular.longitudes.size(); k++)
            for(UInt j=0; j<var.countColumns; j++)
              data.at(idEpoch)(i*griddedDataRectangular.longitudes.size()+k, idxCol+j) = values(var.incLat*i+var.incLon*k+var.incCol*j);
        idxCol += var.countColumns;
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
