/***********************************************/
/**
* @file netCdfInfo.cpp
*
* @brief Content information of a NetCDF file.
*
* @author Torsten Mayer-Guerr
* @date 2020-09-03
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Print content information of a NetCDF file like
dimensions, variables and attributes.

See also \program{NetCdf2GridRectangular}, \program{GridRectangular2NetCdf}.
)";

/***********************************************/

#include "programs/program.h"
#include "inputOutput/fileNetCdf.h"

/***** CLASS ***********************************/

/** @brief Content information of a NetCDF file.
* @ingroup programsConversionGroup */
class NetCdfInfo
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(NetCdfInfo, SINGLEPROCESS, "Content information of a NetCDF file", Conversion)

/***********************************************/

void NetCdfInfo::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameIn;

    readConfig(config, "inputfileNetCdf", fileNameIn, Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

#ifdef NOLIB_NETCDF
    throw(Exception("Compiled without NetCDF library"));
#else
    // open netCDF file
    // ----------------
    logStatus<<"read netCDF file <"<<fileNameIn<<">"<<Log::endl;
    NetCdf::InFile file(fileNameIn);

    logInfo<<"  global attributes:"<<Log::endl;
    auto attributes = file.attributes();
    for(auto &attr : attributes)
      logInfo<<"    - "<<attr.name()<<" = "<<attr.value()<<Log::endl;

    logInfo<<"  dimensions:"<<Log::endl;
    auto dimensions = file.dimensions();
    for(auto &dim : dimensions)
      logInfo<<"    - "<<dim.name()<<" = "<<dim.length()<<Log::endl;

    auto variables = file.variables();
    for(auto &var : variables)
    {
      std::stringstream ss;
      auto dimensions = var.dimensions();
      if(dimensions.size())
      {
        ss<<"  variable: "<<var.name()<<"(";
        ss<<dimensions.at(0).name();
        for(UInt i=1; i<dimensions.size(); i++)
          ss<<", "<<dimensions.at(i).name();
        ss<<")";
      }
      else
        ss<<"  - "<<var.name();
      logInfo<<ss.str()<<Log::endl;

      auto attributes = var.attributes();
      for(auto &attr : attributes)
        logInfo<<"    - "<<attr.name()<<" value = "<<attr.value()<<Log::endl;
    }
#endif
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
