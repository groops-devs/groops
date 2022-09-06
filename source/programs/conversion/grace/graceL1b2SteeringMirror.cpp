/***********************************************/
/**
* @file graceL1b2SteeringMirror.cpp
*
* @brief Read GRACE L1B data.
*
* @author Andreas Kvas
* @date 2019-11-19
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program converts GRACE-FO Steering Mirror output (LSM1B) to an \file{instrument file (STARCAMERA)}{instrument}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "fileGrace.h"

/***** CLASS ***********************************/

/** @brief Read GRACE L1B data.
* @ingroup programsConversionGroup */
class GraceL1b2SteeringMirror
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GraceL1b2SteeringMirror, SINGLEPROCESS, "read GRACE L1B data (LSM1B)", Conversion, Grace, Instrument)

/***********************************************/

void GraceL1b2SteeringMirror::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName              fileNameOut;
    std::vector<FileName> fileNameIn;

    readConfig(config, "outputfileStarCamera", fileNameOut, Config::OPTIONAL, "", "");
    readConfig(config, "inputfile",            fileNameIn,  Config::MUSTSET,  "", "LSM1B");
    if(isCreateSchema(config)) return;

    // =============================================

    logStatus<<"read input files"<<Log::endl;
    StarCameraArc arc;
    for(UInt idFile=0; idFile<fileNameIn.size(); idFile++)
    {
      logStatus<<"read file <"<<fileNameIn.at(idFile)<<">"<<Log::endl;
      UInt numberOfRecords;
      FileInGrace file(fileNameIn.at(idFile), numberOfRecords);

      for(UInt idEpoch=0; idEpoch<numberOfRecords; idEpoch++)
      {
        Int32    seconds;
        Double   frac_seconds;
        Char     gracefo_ID;
        Double   pitch, yaw;
        Byte     qualflg;

        try
        {
          file>>seconds>>frac_seconds>>gracefo_ID>>yaw>>pitch>>FileInGrace::flag(qualflg);
        }
        catch(std::exception &/*e*/)
        {
          // GRACE-FO number of records issue
          logWarning<<arc.back().time.dateTimeStr()<<": file ended at "<<idEpoch<<" of "<<numberOfRecords<<" expected records"<<Log::endl;
          break;
        }

        const Time time = mjd2time(51544.5) + seconds2time(seconds) + seconds2time(frac_seconds*1e-9);
        if(arc.size() && (time <= arc.back().time))
          logWarning<<"epoch("<<time.dateTimeStr()<<") <= last epoch("<<arc.back().time.dateTimeStr()<<")"<<Log::endl;

        StarCameraEpoch epoch;
        epoch.time = time;
        epoch.rotary = rotaryZ(Angle(yaw*1e-6))*rotaryY(Angle(pitch*1e-6));
        arc.push_back(epoch);
      } // for(idEpoch)
    } // for(idFile)

    // =============================================

    logStatus<<"sort epochs"<<Log::endl;
    arc.sort();

    logStatus<<"eliminate duplicates"<<Log::endl;
    const UInt oldSize = arc.size();
    arc.removeDuplicateEpochs(TRUE/*keepFirst*/);
    if(arc.size() < oldSize)
      logInfo<<" "<<oldSize-arc.size()<<" duplicates removed!"<<Log::endl;

    Arc::printStatistics(arc);
    if(arc.size() == 0)
      return;

    if(!fileNameOut.empty())
    {
      logInfo<<"write data to <"<<fileNameOut<<">"<<Log::endl;
      InstrumentFile::write(fileNameOut, arc);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
