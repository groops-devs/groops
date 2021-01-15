/***********************************************/
/**
* @file grailCdr2StarCamera.cpp
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
This program converts orientation data measured by a star camera (SRF to CRF)
from the GRAIL SDS format into \file{instrument file (STARCAMERA)}{instrument}.
For further information see \program{GraceL1b2Accelerometer}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief read CDR GRAIL data.
* @ingroup programsConversionGroup */
class GrailCdr2StarCamera
{
  void readFile(const FileName fileName, StarCameraArc &arc);

public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GrailCdr2StarCamera, SINGLEPROCESS, "read CDR GRAIL data", Conversion, Instrument)

/***********************************************/

void GrailCdr2StarCamera::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName outName;
    std::vector<FileName> fileNames;

    readConfig(config, "outputfileStarCamera", outName,   Config::OPTIONAL, "", "");
    readConfig(config, "inputfile",            fileNames, Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    logStatus<<"read input files"<<Log::endl;
    StarCameraArc arc;
    for(auto &fileName : fileNames)
    {
      logStatus<<"read file <"<<fileName<<">"<<Log::endl;
      readFile(fileName, arc);
    }

    logStatus<<"write data to <"<<outName.str()<<">"<<Log::endl;
    InstrumentFile::write(outName, arc);
    Arc::printStatistics(arc);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GrailCdr2StarCamera::readFile(const FileName fileName, StarCameraArc &arc)
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

    // Header einlesen
    UInt numberOfRecords = 0;
    std::string line;
    for(;;)
    {
      getline(file, line);
      if(line.find("NUMBER OF DATA RECORDS")==0)
        numberOfRecords = String::toInt(line.substr(31, 10));
      // Header fertig ?
      if(line.find("END OF HEADER")==0)
        break;
    }

    // Eigentliche Daten einlesen
    for(UInt i=0; i<numberOfRecords; i++)
    {
      std::string line;
      try
      {
        getline(file, line);
      }
      catch(std::exception &/*e*/)
      {
        //logWarning<<std::endl<<e.what()<<" continue..."<<Log::endl;
        break;
      }

      std::stringstream ss(line);

      // time (TDB) >> receiver time, seconds past 12:00:00 noon 01-Jan-2000
      Int32 seconds;
      ss>>seconds;
      StarCameraEpoch epoch;
      epoch.time = mjd2time(51544.5) + seconds2time(seconds);

      // satellite ID, SCA identification number
      Byte AB;
      Int  SCA;
      ss>>AB>>SCA;

      // quaternions + formal error of quaternions
      Vector q(4);
      Double q_res;
      ss>>q(0)>>q(1)>>q(2)>>q(3)>>q_res;
      epoch.rotary = Rotary3d(q);

      // data quality flags >> not defined yet
      //Byte flag;
      //ss>>flag;

      arc.push_back(epoch);
    }
  }

  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
