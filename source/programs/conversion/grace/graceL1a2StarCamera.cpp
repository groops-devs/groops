/***********************************************/
/**
* @file graceL1a2StarCamera.cpp
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
from the GRACE Level-1A format (SCA1A) to the GROOPS instrument file format.
For further information see \program{GraceL1a2Accelerometer}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "fileGrace.h"

/***** CLASS ***********************************/

/** @brief Read GRACE L1A data.
* @ingroup programsConversionGroup */
class GraceL1a2StarCamera
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GraceL1a2StarCamera, SINGLEPROCESS, "read GRACE L1A data (SCA1A)", Conversion, Grace, Instrument)

/***********************************************/

void GraceL1a2StarCamera::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameSca1, fileNameSca2;
    std::vector<FileName> fileNameIn;

    readConfig(config, "outputfileStarCamera1", fileNameSca1, Config::MUSTSET,  "", "STARCAMERA1A, head 1");
    readConfig(config, "outputfileStarCamera2", fileNameSca2, Config::MUSTSET,  "", "STARCAMERA1A, head 2");
    readConfig(config, "inputfile",             fileNameIn,   Config::MUSTSET,  "", "SCA1A, !GRACE-FO is not working!");
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
        Int32             seconds, microSeconds=0;
        Char              GRACE_id;
        FileInGrace::Int8 sca_id;
        Char              scaDesign;
        Double            q0, q1, q2, q3;
        FileInGrace::Int8 nLocks, nStars;
        FileInGrace::Int8 scaConfig1, scaConfig2, scaConfig3;
        Byte              scaMode;
        Byte              qualityFlag;

        try
        {
          file>>seconds>>/*microSeconds>>*/GRACE_id>>sca_id>>scaDesign>>q0>>q1>>q2>>q3; // microSeconds in GRACE-FO data only
          file>>nLocks>>nStars>>scaConfig1>>scaConfig2>>scaConfig3>>FileInGrace::flag(scaMode)>>FileInGrace::flag(qualityFlag);
        }
        catch(std::exception &/*e*/)
        {
          // GRACE-FO number of records issue
          logWarning<<arc1.back().time.dateTimeStr()<<": file ended at "<<idEpoch<<" of "<<numberOfRecords<<" expected records"<<Log::endl;
          break;
        }

       {
          StarCamera1AEpoch epoch;
          epoch.time = mjd2time(51544.5) + seconds2time(seconds) + seconds2time(1e-6*microSeconds);
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

          if((UInt(sca_id) == 1) && !(qualityFlag & (1 << 0)) && (UInt(scaConfig1) < 7)) // invalid data is removed
            arc1.push_back(epoch);

          if((UInt(sca_id) == 2) && !(qualityFlag & (1 << 0)) && (UInt(scaConfig1) < 7)) // invalid data is removed
            arc2.push_back(epoch);
        }
      } // for(idEpoch)
    } // for(idFile)

    // =============================================

    if(!fileNameSca1.empty())
    {
      logStatus<<"Star camera head 1:"<<Log::endl;
      Arc::printStatistics(arc1);
      if(arc1.size())
      {
        logStatus<<"write data to <"<<fileNameSca1<<">"<<Log::endl;
        InstrumentFile::write(fileNameSca1, arc1);
      }
    }

    if(!fileNameSca2.empty())
    {
      logStatus<<"Star camera head 2:"<<Log::endl;
      Arc::printStatistics(arc2);
      if(arc2.size())
      {
        logStatus<<"write data to <"<<fileNameSca2<<">"<<Log::endl;
        InstrumentFile::write(fileNameSca2, arc2);
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
