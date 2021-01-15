/***********************************************/
/**
* @file jplAscii2Ephemerides.cpp
*
* @brief Read JPL DExxx ephemerides.
*
* @author Torsten Mayer-Guerr
* @date 2020-08-11
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Read JPL DExxx (ASCII) ephemerides.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "inputOutput/file.h"
#include "files/fileEphemerides.h"

/***** CLASS ***********************************/

/** @brief Read JPL DExxx ephemerides.
* @ingroup programsConversionGroup */
class JplAscii2Ephemerides
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(JplAscii2Ephemerides, SINGLEPROCESS, "JPL DExxx (ASCII) ephemerides", Conversion)

/***********************************************/

void JplAscii2Ephemerides::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName              fileNameOut, fileNameHeader;
    std::vector<FileName> fileNamesIn;

    readConfig(config, "outputfileEphemerides", fileNameOut,    Config::MUSTSET,  "JPL_DE432.dat",  "");
    readConfig(config, "inputfileHeader",       fileNameHeader, Config::MUSTSET,  "header.432_571", "");
    readConfig(config, "inputfileData",         fileNamesIn,    Config::MUSTSET,  "ascp01950.432",  "");
    if(isCreateSchema(config)) return;

    // ==============================

    auto readDouble = [](InFile &file)
    {
      std::string str;
      file>>str;
      return String::toDouble(str);
    };

    logStatus<<"read header file <"<<fileNameHeader<<">"<<Log::endl;
    Double            earthMoonRatio = 0;
    std::vector<UInt> index(15), subintervals(15), components(15), degree(15);

    InFile file(fileNameHeader);
    std::string line;
    while(std::getline(file, line))
    {
      if(String::startsWith(line, "GROUP   1041")) // constants
      {
        UInt count;
        file>>count;
        for(UInt i=0; i<11; i++)
          earthMoonRatio = readDouble(file);
      }

      if(String::startsWith(line, "GROUP   1050")) // index table
      {
        for(UInt i=0; i<index.size(); i++)
          file>>index.at(i);
        for(UInt i=0; i<degree.size(); i++)
          file>>degree.at(i);
        for(UInt i=0; i<subintervals.size(); i++)
          file>>subintervals.at(i);
        components = {3,  // MERCURY
                      3,  // VENUS
                      3,  // EARTH
                      3,  // MARS
                      3,  // JUPITER
                      3,  // SATURN
                      3,  // URANUS
                      3,  // NEPTUNE
                      3,  // PLUTO
                      3,  // MOON
                      3,  // SUN
                      2,  // NUTATION
                      3,  // LIBRATION
                      3,  // Lunar mantle angular vel
                      1}; // TT - TDB
        for(UInt i=0; i<degree.size(); i++)
          if(degree.at(i))
            degree.at(i)--;
        for(UInt i=0; i<subintervals.size(); i++)
          if(subintervals.at(i) == 0)
            components.at(i) = 0;
      }
    }

    // ==============================

    std::vector<Time> times;
    std::vector<std::vector<std::vector<Matrix>>> coefficients;
    for(FileName &filenameIn : fileNamesIn)
    {
      try
      {
        logStatus<<"read file <"<<filenameIn<<">"<<Log::endl;
        InFile file(filenameIn);
        for(;;)
        {
          UInt idx, count;
          file>>idx>>count;
          if(file.eof())
            break;

          Time timeStart = timeTT2GPS(mjd2time(readDouble(file)-2400000.5));
          Time timeEnd   = timeTT2GPS(mjd2time(readDouble(file)-2400000.5));
          logInfo<<"["<<timeStart.dateTimeStr()<<", "<<timeEnd.dateTimeStr()<<") days = "<<(timeEnd-timeStart).mjd()<<Log::endl;

          if(!times.size())
            times.push_back(timeStart);

          UInt idInterval = std::distance(times.begin(), std::find(times.begin(), times.end(), timeStart));
          if(idInterval+1 < times.size()) // interval already exists?
          {
            if(times.at(idInterval+1) != timeEnd)
              throw(Exception("expected "+times.at(idInterval+1).dateTimeStr()+" but get "+timeEnd.dateTimeStr()));
          }
          else
          {
            if(idInterval+1 != times.size())
              throw(Exception("expected "+times.back().dateTimeStr()+" but get "+timeStart.dateTimeStr()));
            times.push_back(timeEnd);
            coefficients.resize(times.size()-1);
            coefficients.back().resize(subintervals.size());
          }

          for(UInt idBody=0; idBody<subintervals.size(); idBody++)
          {
            coefficients.at(idInterval).at(idBody).resize(subintervals.at(idBody), Matrix(components.at(idBody), degree.at(idBody)+1));
            for(UInt idSub=0; idSub<subintervals.at(idBody); idSub++)
              for(UInt i=0; i<components.at(idBody); i++)
                for(UInt n=0; n<=degree.at(idBody); n++)
                  coefficients.at(idInterval).at(idBody).at(idSub)(i, n) = readDouble(file); // km -> m
          }
        }
      }
      catch(std::exception &e)
      {
        logWarning<<std::endl<<e.what()<<" continue..."<<Log::endl;
        continue;
      }
    } // for(inputFiles)

    // ==============================

    // km -> m
    // -------
    for(UInt idInterval=0; idInterval<coefficients.size(); idInterval++)
      for(UInt idBody=0; idBody<std::min(std::size_t(11), subintervals.size()); idBody++)
        for(UInt idSub=0; idSub<subintervals.at(idBody); idSub++)
          coefficients.at(idInterval).at(idBody).at(idSub) *= 1000;

    // ==============================

    // write results
    // -------------
    logStatus<<"write ephemerides to file <"<<fileNameOut<<">"<<Log::endl;
    writeFileEphemerides(fileNameOut, earthMoonRatio, times, subintervals, components, degree, coefficients);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/
