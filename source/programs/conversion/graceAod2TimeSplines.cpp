/***********************************************/
/**
* @file graceAod2TimeSplines.cpp
*
* @brief Convert AOD1B dealiasing data into linear splines.
*
* @author Torsten Mayer-Guerr
* @date 2005-01-19
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program converts the atmospheric and ocean de-aliasing product (AOD1B)
from the GRACE SDS format into \file{time spline files}{timeSplinesGravityField}.
Multiple \config{inputfile}s must be given in the correct time order.
A linear method is assumed for the interpolation between the given points in time.

The GRACE SDS format is described in "AOD1B Product Description Document"
given at \url{http://podaac.jpl.nasa.gov/grace/documentation.html}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "inputOutput/file.h"
#include "files/fileTimeSplinesGravityfield.h"

/***** CLASS ***********************************/

/** @brief Convert AOD1B dealiasing data into linear splines.
* @ingroup programsConversionGroup */
class GraceAod2TimeSplines
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GraceAod2TimeSplines, SINGLEPROCESS, "convert AOD1B dealiasing data into linear splines", Conversion, TimeSplines)
GROOPS_RENAMED_PROGRAM(GraceCsrAod2TimeSplines, GraceAod2TimeSplines, date2time(2020, 6, 14))

/***********************************************/

void GraceAod2TimeSplines::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName outputName, atmosName, oceanName, obaName, miscName;
    std::vector<FileName> fileNameIn;

    readConfig(config, "outputfileDealiasing",     outputName, Config::OPTIONAL, "", "");
    readConfig(config, "outputfileAtmosphere",     atmosName,  Config::OPTIONAL, "", "");
    readConfig(config, "outputfileOcean",          oceanName,  Config::OPTIONAL, "", "");
    readConfig(config, "outputfileBottomPressure", obaName,    Config::OPTIONAL, "", "");
    readConfig(config, "outputfileMisc",           miscName,   Config::OPTIONAL, "", "");
    readConfig(config, "inputfile",                fileNameIn, Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    Double GM = DEFAULT_GM;
    Double R  = DEFAULT_R;
    std::vector<Time>   timeList, timeAtmosList, timeOceanList, timeObaList, timeMiscList;
    std::vector<Matrix> cnmList,      snmList;
    std::vector<Matrix> cnmAtmosList, snmAtmosList;
    std::vector<Matrix> cnmOceanList, snmOceanList;
    std::vector<Matrix> cnmObaList,   snmObaList;
    std::vector<Matrix> cnmMiscList,  snmMiscList;

    for(UInt idFile=0; idFile<fileNameIn.size(); idFile++)
    {
      try
      {
        logStatus<<"read file <"<<fileNameIn.at(idFile)<<">"<<Log::endl;
        InFile file(fileNameIn.at(idFile));

        // Header
        std::string line;
        UInt dataCount = 0;
        UInt version   = 9999;
        UInt degree    = 100;
        while(!file.eof())
        {
          std::getline(file, line);
          if(line.find("CONSTANT GM")==0)
            GM = String::toDouble(line.substr(31, 22));
          if(line.find("CONSTANT A")==0)
            R = String::toDouble(line.substr(31, 22));
          if(line.find("MAXIMUM DEGREE")==0)
            degree = String::toInt(line.substr(32, 3));
          if(line.find("NUMBER OF DATA SETS")==0)
            dataCount = String::toInt(line.substr(31, 22));
          if(line.find("SOFTWARE VERSION")==0)
            version = String::toInt(line.substr(50, 2));
          if(line.find("END OF HEADER")==0)
            break;
        }

        if(version==9999)
          logWarning<<"No SOFTWARE VERSION header record found"<<Log::endl;
        if(dataCount==0)
          logWarning<<"No NUMBER OF DATA SETS header record found"<<Log::endl;

        for(UInt k=0; k<dataCount; k++)
        {
          // Data Header
          std::getline(file, line);

          Int year, month, day, hour;
          if(version==0)
          {
            year  = String::toInt(line.substr(39, 4));
            month = String::toInt(line.substr(44, 2));
            day   = String::toInt(line.substr(47, 2));
            hour  = String::toInt(line.substr(50, 2));
          }
          else
          {
            year  = String::toInt(line.substr(37, 4));
            month = String::toInt(line.substr(42, 2));
            day   = String::toInt(line.substr(45, 2));
            hour  = String::toInt(line.substr(48, 2));
          }
          Time time = date2time(year, month, day, hour, 0, 0);

          std::string type;
          if(version==0)
            type = "glo";
          else
            type = line.substr(65,3);

          Matrix cnm(degree+1, Matrix::TRIANGULAR, Matrix::LOWER);
          Matrix snm(degree+1, Matrix::TRIANGULAR, Matrix::LOWER);

          if(type != "tst")
          {
            for(UInt n=0; n<=degree; n++)
              for(UInt m=0; m<=n; m++)
              {
                UInt n2,m2;
                file>>n2>>m2;
                file>>cnm(n2,m2)>>snm(n2,m2);
                std::getline(file, line); // get rest of line
              }
          }

          if(type == "glo")
          {
            timeList.push_back(time);
            cnmList.push_back(cnm);
            snmList.push_back(snm);
          }
          else if(type == "atm")
          {
            timeAtmosList.push_back(time);
            cnmAtmosList.push_back(cnm);
            snmAtmosList.push_back(snm);
          }
          else if(type == "ocn")
          {
            timeOceanList.push_back(time);
            cnmOceanList.push_back(cnm);
            snmOceanList.push_back(snm);
          }
          else if(type == "oba")
          {
            timeObaList.push_back(time);
            cnmObaList.push_back(cnm);
            snmObaList.push_back(snm);
          }
          else
          {
            logInfo<<"Type: "<<type<<", "<<time.dateTimeStr()<<Log::endl;
            timeMiscList.push_back(time);
            cnmMiscList.push_back(cnm);
            snmMiscList.push_back(snm);
          }
        }
      }
      catch(std::exception &e)
      {
        logError<<e.what()<<": continue..."<<Log::endl;
      }
    }

    if(!isRegular(timeList) || !isRegular(timeAtmosList) || !isRegular(timeOceanList) || !isRegular(timeObaList) || !isRegular(timeMiscList))
    {
      const Time sampling = medianSampling(timeList);
      for(UInt i=0; i<timeList.size()-1; i++)
        if((timeList.at(i+1)-timeList.at(i)) > sampling)
          logWarning<<"gap between "<<timeList.at(i).dateTimeStr()<<" and "<<timeList.at(i+1).dateTimeStr()<<Log::endl;
      throw(Exception("Spline time series is not regular."));
    }

    // write data
    // ----------
    if(!outputName.empty())
    {
      logStatus<<"write dealiasing data to <"<<outputName.str()<<">"<<Log::endl;
      writeFileTimeSplinesGravityfield(outputName, GM, R, 1/*splineDegree*/, timeList, cnmList, snmList);
    }

    if(!atmosName.empty())
    {
      logStatus<<"write atmosphere data to <"<<atmosName.str()<<">"<<Log::endl;
      writeFileTimeSplinesGravityfield(atmosName, GM, R, 1/*splineDegree*/, timeAtmosList, cnmAtmosList, snmAtmosList);
    }

    if(!oceanName.empty())
    {
      logStatus<<"write ocean data to <"<<oceanName.str()<<">"<<Log::endl;
      writeFileTimeSplinesGravityfield(oceanName, GM, R, 1/*splineDegree*/, timeOceanList, cnmOceanList, snmOceanList);
    }

    if(!obaName.empty())
    {
      logStatus<<"write bottom pressure data to <"<<obaName.str()<<">"<<Log::endl;
      writeFileTimeSplinesGravityfield(obaName, GM, R, 1/*splineDegree*/, timeObaList, cnmObaList, snmObaList);
    }

    if(!miscName.empty())
    {
      logStatus<<"write misc data to <"<<miscName.str()<<">"<<Log::endl;
      writeFileTimeSplinesGravityfield(miscName, GM, R, 1/*splineDegree*/, timeMiscList, cnmMiscList, snmMiscList);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
