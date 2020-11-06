/***********************************************/
/**
* @file iersC04IAU2000EarthOrientationParameter.cpp
*
* @brief Read Earth Orientation Parameter.
*
* @author Torsten Mayer-Guerr
* @date 2005-01-14
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Read a IERS Earth orientation data C04 (IAU2000A) file
and write it as \configFile{outputfileEOP}{earthOrientationParameter}.
)";

/***********************************************/

#include "programs/program.h"
#include "inputOutput/file.h"
#include "files/fileEarthOrientationParameter.h"

/***** CLASS ***********************************/

/** @brief Read Earth Orientation Parameter.
* @ingroup programsConversionGroup */
class IersC04IAU2000EarthOrientationParameter
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(IersC04IAU2000EarthOrientationParameter, SINGLEPROCESS, "read Earth Orientation Parameter", Conversion)
GROOPS_RENAMED_PROGRAM(Eop2003file, IersC04IAU2000EarthOrientationParameter, date2time(2020, 9, 8))

/***********************************************/

void IersC04IAU2000EarthOrientationParameter::run(Config &config)
{
  try
  {
    FileName inName, outName;
    Time timeStart = Time();
    Time timeEnd   = date2time(9999,1,1,0,0,0);

    readConfig(config, "outputfileEOP", outName,   Config::MUSTSET,  "", "");
    readConfig(config, "inputfile",     inName,    Config::MUSTSET,  "", "");
    readConfig(config, "timeStart",     timeStart, Config::OPTIONAL, "", "");
    readConfig(config, "timeEnd",       timeEnd,   Config::OPTIONAL, "", "");
    if(isCreateSchema(config)) return;

    // Erdrotationsparameter einlesen
    // ------------------------------
    logStatus<<"read input file"<<Log::endl;
    std::vector<Double> mjd, xp, yp, du, ld, dX, dY;

    InFile file(inName);

    // Header
    std::string line;
    for(UInt i=0; i<15; i++)
      std::getline(file, line);

    while(std::getline(file, line))
    {
      std::stringstream ss(line);
      ss.exceptions(std::ios::badbit | std::ios::failbit);

      UInt   day, month, year;
      Double _mjd, _xp, _yp, _du, _ld, _dX, _dY;
      ss>>year>>month>>day>>_mjd>>_xp>>_yp>>_du>>_ld>>_dX>>_dY;

      const Time time = mjd2time(_mjd);
      if((time < timeStart) || (time > timeEnd))
        continue;

      mjd.push_back(_mjd);
      xp.push_back(_xp);
      yp.push_back(_yp);
      du.push_back(_du);
      ld.push_back(_ld);
      dX.push_back(_dX);
      dY.push_back(_dY);
    }

    logInfo<<"  count = "<<mjd.size()<<Log::endl;
    logInfo<<"  start = "<<mjd2time(mjd.at(0)).dateTimeStr()<<Log::endl;
    logInfo<<"  end   = "<<mjd2time(mjd.back()).dateTimeStr()<<Log::endl;

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
    logStatus<<"write EOPs <"<<outName<<">"<<Log::endl;
    writeFileEarthOrientationParameter(outName, EOP);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
