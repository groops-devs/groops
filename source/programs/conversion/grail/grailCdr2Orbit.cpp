/***********************************************/
/**
* @file grailCdr2Orbit.cpp
*
* @brief read CDR GRAIL data.
*
* @author Beate Klinger
* @date 2013-01-22
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program converts the orbit from the GRAIL SDS format into  \file{instrument file (ORBIT)}{instrument}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief read CDR GRAIL data.
* @ingroup programsConversionGroup */
class GrailCdr2Orbit
{
  void readFile(const FileName fileName, OrbitArc &arc);

public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GrailCdr2Orbit, SINGLEPROCESS, "read CDR GRAIL data", Conversion, Orbit, Instrument)

/***********************************************/

void GrailCdr2Orbit::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName outName;
    std::vector<FileName> fileNames;

    readConfig(config, "outputfileOrbit", outName,   Config::MUSTSET,  "", "");
    readConfig(config, "inputfile",       fileNames, Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    logStatus<<"read input files"<<Log::endl;
    OrbitArc orbit;
    for(auto &fileName : fileNames)
    {
      logStatus<<"read file <"<<fileName<<">"<<Log::endl;
      readFile(fileName, orbit);
    }

    logStatus<<"write data to <"<<outName<<">"<<Log::endl;
    InstrumentFile::write(outName, orbit);
    Arc::printStatistics(orbit);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GrailCdr2Orbit::readFile(const FileName fileName, OrbitArc &arc)
{
  try
  {
    InFile file(fileName);
    if(!file.good())
    {
      logWarning<<"cannot open file: "<<fileName.str()<<", continue..."<<Log::endl;
      return;
    }
    file.exceptions(std::ios::badbit|std::ios::failbit);

    // Header
    UInt numberOfRecords = 0;
    std::string line;
    for(;;)
    {
      getline(file, line);
      if(line.find("NUMBER OF DATA RECORDS")==0)
        numberOfRecords = String::toInt(line.substr(31, 10));
      if(line.find("END OF HEADER")==0)
        break;
    }

    // Data
    for(UInt i=0; i<numberOfRecords; i++)
    {
      std::string line;
      try
      {
        std::getline(file, line);
      }
      catch(std::exception &/*e*/)
      {
        //logWarning<<std::endl<<e.what()<<" continue..."<<Log::endl;
        break;
      }
      std::stringstream ss(line);

      // time (TDB) >> seconds past 12:00:00 noon 01-Jan-2000
      UInt seconds;
      ss>>seconds;
      OrbitEpoch epoch;
      epoch.time = timeTT2GPS(mjd2time(51544.5)+seconds2time(seconds));

      // satellite ID, coordinate reference frame
      Int32 AB;
      Byte  IE; //Byte
      ss>>AB>>IE;

      // position + formal error in position
      Double x_error, y_error, z_error;
      ss>>epoch.position.x()>>epoch.position.y()>>epoch.position.z();
      ss>>x_error>>y_error>>z_error;

      // velocity + formalerror in velocity
      Double dx_error, dy_error, dz_error;
      ss>>epoch.velocity.x()>>epoch.velocity.y()>>epoch.velocity.z();
      ss>>dx_error>>dy_error>>dz_error;

      // data quality flags >> data quality flags not defined!
      //Int flag;
      //ss>>flag;
      //if(flag!=0)
      //  continue;

      arc.push_back(epoch);
    }
  }

  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
