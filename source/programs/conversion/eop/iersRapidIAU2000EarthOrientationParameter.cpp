/***********************************************/
/**
* @file iersRapidIAU2000EarthOrientationParameter.cpp
*
* @brief Read Earth Orientation Parameter from rapid file format.
*
* @author Andreas KVas
* @date 2016-03-21
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Read a IERS Earth orientation rapid data and prediction file (IAU2000)
and write it as \configFile{outputfileEOP}{earthOrientationParameter}.
)";

/***********************************************/

#include "programs/program.h"
#include "inputOutput/file.h"
#include "files/fileEarthOrientationParameter.h"

/***** CLASS ***********************************/

/** @brief Read Earth Orientation Parameter.
* @ingroup programsConversionGroup */
class IersRapidIAU2000EarthOrientationParameter
{
  Double getDouble(std::string s, const std::vector<Double>& vals, double f=1.0) const;
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(IersRapidIAU2000EarthOrientationParameter, SINGLEPROCESS, "read Earth Orientation Parameter", Conversion)
GROOPS_RENAMED_PROGRAM(IersRapidEop2Earthrotation, IersRapidIAU2000EarthOrientationParameter, date2time(2020, 9, 8))

/***********************************************/

void IersRapidIAU2000EarthOrientationParameter::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName inName, outName;
    Time timeStart = Time();
    Time timeEnd   = date2time(9999,1,1);

    readConfig(config, "outputfileEOP", outName,   Config::MUSTSET,  "", "");
    readConfig(config, "inputfile",     inName,    Config::MUSTSET,  "", "");
    readConfig(config, "timeStart",     timeStart, Config::OPTIONAL, "", "");
    readConfig(config, "timeEnd",       timeEnd,   Config::OPTIONAL, "", "");
    if(isCreateSchema(config)) return;

    // Erdrotationsparameter einlesen
    // ------------------------------
    logStatus<<"read input file <"<<inName<<">"<<Log::endl;
    InFile file(inName);

    std::vector<Double> mjd, xp, yp, du, ld, dX, dY;
    std::string line;
    while(std::getline(file, line))
    {
      if(line.empty())
        continue;

      Double _mjd, _xp, _yp, _du, _dX, _dY, _ld;
      if(!(std::stringstream(line.substr(7, 8))>>_mjd))
        throw(Exception("Error reading MJD from EOP file."));

      Time time = mjd2time(_mjd);
      if((time<timeStart)||(time>timeEnd))
        continue;

      _xp = getDouble(line.substr( 18,  9),xp);
      _yp = getDouble(line.substr( 37,  9),yp);
      _du = getDouble(line.substr( 58, 10),du);
      _ld = getDouble(line.substr( 79,  7),ld,1e3);
      _dX = getDouble(line.substr( 97,  9),dX,1e3);
      _dY = getDouble(line.substr(116,  9),dY,1e3);

      mjd.push_back(_mjd);
      xp.push_back(_xp);
      yp.push_back(_yp);
      du.push_back(_du);
      ld.push_back(_ld*1e-3);
      dX.push_back(_dX*1e-3);
      dY.push_back(_dY*1e-3);
    }

    logInfo<<"  count = "<<mjd.size()<<Log::endl;
    if(mjd.size())
    {
      logInfo<<"  start = "<<mjd2time(mjd.at(0)).dateTimeStr()<<Log::endl;
      logInfo<<"  end   = "<<mjd2time(mjd.back()).dateTimeStr()<<Log::endl;
    }

    Matrix EOP(mjd.size(), 7);
    for(UInt i=0; i<EOP.rows(); i++)
    {
      EOP(i,0) = mjd.at(i);
      EOP(i,1) = xp.at(i);
      EOP(i,2) = yp.at(i);
      EOP(i,3) = du.at(i);
      EOP(i,4) = ld.at(i);
      EOP(i,5) = dX.at(i);
      EOP(i,6) = dY.at(i);
    }

    // Save to file
    // ------------
    logStatus<<"write EOPs"<<Log::endl;
    writeFileEarthOrientationParameter(outName, EOP);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double IersRapidIAU2000EarthOrientationParameter::getDouble(std::string s, const std::vector<Double>& vals, double f) const
{
  Double val;
  if(std::stringstream(s)>>val)
    return val;
  else if(vals.size())
      return vals.back()*f;
  else
    throw(Exception("Error reading value from EOP file."));
}
