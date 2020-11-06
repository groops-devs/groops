/***********************************************/
/**
* @file graceL1A2StarCamera.cpp
*
* @brief Read GRACE L1A data.
*
* @author Beate Klinger
* @date 2017-01-03
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program converts orientation data measured by the star cameras
from the GRACE Level-1A format to the GROOPS instrument file format.
For further information see \program{GraceL1A2Accelerometer}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "fileGrace.h"

/***** CLASS ***********************************/

/** @brief Read GRACE L1A data.
* @ingroup programsConversionGroup */
class GraceL1A2StarCamera
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(GraceL1A2StarCamera, SINGLEPROCESS, "read GRACE L1A data", Conversion, Grace, Instrument)

/***********************************************/

void GraceL1A2StarCamera::run(Config &config)
{
  try
  {
    FileName fileNameSca1, fileNameSca2;
    std::vector<FileName> fileNameIn;

    readConfig(config, "outputfileStarCamera1", fileNameSca1, Config::MUSTSET,  "", "");
    readConfig(config, "outputfileStarCamera2", fileNameSca2, Config::MUSTSET,  "", "");
    readConfig(config, "inputfile",             fileNameIn,      Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    // =============================================

    logStatus<<"read input files"<<Log::endl;
    Arc arc1, arc2;
    for(UInt idFile=0; idFile<fileNameIn.size(); idFile++)
    {
      logStatus<<"read file <"<<fileNameIn.at(idFile)<<">"<<Log::endl;
      UInt numberOfRecords;
      FileInGrace file(fileNameIn.at(idFile), numberOfRecords);

      for(UInt idEpoch=0; idEpoch<numberOfRecords; idEpoch++)
      {
        Int32  seconds;
        Byte   GRACE_id, sca_id, scaDesign;
        Double q0, q1, q2, q3;
        Byte   nLocks, nStars;
        Byte   scaConfig1, scaConfig2, scaConfig3, scaMode, qualityFlag;

        file>>seconds>>GRACE_id>>sca_id>>scaDesign>>q0>>q1>>q2>>q3>>nLocks>>nStars>>scaConfig1>>scaConfig2>>scaConfig3>>scaMode>>FileInGrace::flag(qualityFlag);

        const Time time = mjd2time(51544.5) + seconds2time(seconds);
        if(arc1.size() && (UInt(sca_id) == 1) && (time <= arc1.at(arc1.size()-1).time))
          logWarning<<"SCA1: epoch("<<time.dateTimeStr()<<") <= last epoch("<<arc1.at(arc1.size()-1).time.dateTimeStr()<<")"<<Log::endl;
        if(arc2.size() && (UInt(sca_id) == 2) && (time <= arc2.at(arc2.size()-1).time))
          logWarning<<"SCA2: epoch("<<time.dateTimeStr()<<") <= last epoch("<<arc2.at(arc2.size()-1).time.dateTimeStr()<<")"<<Log::endl;

        {
          StarCamera1AEpoch epoch;
          epoch.time = mjd2time(51544.5) + seconds2time(seconds);
          epoch.rcvTime = seconds;
          epoch.epsTime = 0.0;
          if(scaDesign == 'P') // primary = 1, secondary = 2
            epoch.scaDesign = 1;
          if(scaDesign == 'S')
            epoch.scaDesign = 2;
          epoch.q0      = q0;
          epoch.q1      = q1;
          epoch.q2      = q2;
          epoch.q3      = q3;
          epoch.nLocks  = UInt(nLocks);
          epoch.nStars  = UInt(nStars);

          if((UInt(sca_id) == 1) && (Bool(qualityFlag & (1 << 0)) == 0) && (UInt(scaConfig1) < 7)) // invalid data is removed
          arc1.push_back(epoch);

          if((UInt(sca_id) == 2) && (Bool(qualityFlag & (1 << 0)) == 0) && (UInt(scaConfig1) < 7)) // invalid data is removed
          arc2.push_back(epoch);
        }
      } // for(idEpoch)
    } // for(idFile)

    // =============================================

    logInfo<<"Star camera head 1:"<<Log::endl;
    Arc::printStatistics(arc1);

    // remove duplicates
    arc1.sort();
    UInt countDuplicates1 = 0;
    for(UInt idEpoch=1; idEpoch<arc1.size(); idEpoch++)
    {
      if(arc1.at(idEpoch).time == arc1.at(idEpoch-1).time)
      {
        arc1.remove(idEpoch--);
        countDuplicates1++;
      }
    }
    logInfo<<"  duplicates:      "<<countDuplicates1<<Log::endl;

    if(!fileNameSca1.empty())
    {
      logInfo<<"write data to <"<<fileNameSca1<<">"<<Log::endl;
      InstrumentFile::write(fileNameSca1, arc1);
    }

    // =============================================

    logInfo<<"Star camera head 2:"<<Log::endl;
    Arc::printStatistics(arc2);

    // remove duplicates
    arc2.sort();
    UInt countDuplicates2 = 0;
    for(UInt idEpoch=1; idEpoch<arc2.size(); idEpoch++)
    {
      if(arc2.at(idEpoch).time == arc2.at(idEpoch-1).time)
      {
        arc2.remove(idEpoch--);
        countDuplicates2++;
      }
    }
    logInfo<<"  duplicates:      "<<countDuplicates2<<Log::endl;

    if(!fileNameSca2.empty())
    {
      logInfo<<"write data to <"<<fileNameSca2<<">"<<Log::endl;
      InstrumentFile::write(fileNameSca2, arc2);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
