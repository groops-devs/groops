/***********************************************/
/**
* @file igs2EarthOrientationParameter.cpp
*
* @brief Earth Orientation Parameter from IGS daily file.
*
* @author Andreas Kvas
* @date 2015-04-20
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Read Rapid Earth Orientation Parameter from IGS daily file
and write it as \configFile{outputfileEOP}{earthOrientationParameter}.
)";

/***********************************************/

#include "programs/program.h"
#include "inputOutput/file.h"
#include "files/fileEarthOrientationParameter.h"

/***** CLASS ***********************************/

/** @brief Earth Orientation Parameter from IGS daily file.
* @ingroup programsConversionGroup */
class Igs2EarthOrientationParameter
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Igs2EarthOrientationParameter, SINGLEPROCESS, "Earth Orientation Parameter from IGS daily file", Conversion)
GROOPS_RENAMED_PROGRAM(IGS2EarthRotation, Igs2EarthOrientationParameter, date2time(2020, 9, 8))

/***********************************************/

void Igs2EarthOrientationParameter::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName              outName;
    std::vector<FileName> fileNameIn;
    Time timeStart, timeEnd = date2time(9999,1,1);

    readConfig(config, "outputfileEOP", outName,    Config::MUSTSET,  "", "");
    readConfig(config, "inputfile",     fileNameIn, Config::MUSTSET,  "", "");
    readConfig(config, "timeStart",     timeStart,  Config::OPTIONAL, "", "");
    readConfig(config, "timeEnd",       timeEnd,    Config::OPTIONAL, "", "");
    if(isCreateSchema(config)) return;

    logStatus<<"read input files"<<Log::endl;
    std::vector<Double> mjd, xp, yp, du, ld, dX, dY;
    for(UInt idFile=0; idFile<fileNameIn.size(); idFile++)
    {
      logStatus<<"read file <"<<fileNameIn.at(idFile)<<">"<<Log::endl;
      InFile file(fileNameIn.at(idFile));

      std::string line;
      for(UInt i=0; i<9; i++)
        std::getline(file, line); // header

      Double _mjd, _xp, _yp, _du, _ld;
      file>>_mjd>>_xp>>_yp>>_du>>_ld;

      if((_mjd < timeStart.mjd()) || (_mjd > timeEnd.mjd()))
        continue;

      mjd.push_back(_mjd);
      xp.push_back(_xp*1e-6);
      yp.push_back(_yp*1e-6);
      du.push_back(_du*1e-7);
      ld.push_back(_ld*1e-7);
      dX.push_back(0.0);
      dY.push_back(0.0);
    } // for(idFile)


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

    logStatus<<"write EOPs"<<Log::endl;
    writeFileEarthOrientationParameter(outName, EOP);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
