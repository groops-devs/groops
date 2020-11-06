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
This program converts GRACE-FO Steering Mirror output to an \file{instrument file (STARCAMERA)}{instrument}.
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
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(GraceL1b2SteeringMirror, SINGLEPROCESS, "read GRACE L1B data", Conversion, Grace, Instrument)

/***********************************************/

void GraceL1b2SteeringMirror::run(Config &config)
{
  try
  {
    FileName              fileNameOut;
    std::vector<FileName> fileNameIn;

    readConfig(config, "outputfileStarCamera",        fileNameOut,       Config::OPTIONAL, "", "");
    readConfig(config, "inputfile",                   fileNameIn,        Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    // =============================================

    logStatus<<"read input files"<<Log::endl;
    StarCameraArc arcPitchYaw;
    for(UInt idFile=0; idFile<fileNameIn.size(); idFile++)
    {
      logStatus<<"read file <"<<fileNameIn.at(idFile)<<">"<<Log::endl;
      UInt numberOfRecords;
      FileInGrace file(fileNameIn.at(idFile), numberOfRecords);

      for(UInt idEpoch=0; idEpoch<numberOfRecords; idEpoch++)
      {
        Int32    seconds;
        Double   frac_seconds, pitch, yaw;
        Char     gracefo_ID;
        Byte     qualflg;

        try //This block is added for GRACE-FO number of records issue
        {
          file>>seconds>>frac_seconds;
        }
        catch(std::exception &/*e*/)
        {
          logWarning<<arcPitchYaw.at(arcPitchYaw.size()-1).time.dateTimeStr()<<": file ended at "<<idEpoch<<" of "<<numberOfRecords<<" expected records"<<Log::endl;
          break;
        }

        file>>gracefo_ID>>yaw>>pitch;
        file>>FileInGrace::flag(qualflg);

        const Time time = mjd2time(51544.5) + seconds2time(seconds) + seconds2time(frac_seconds*1e-9);
        if(arcPitchYaw.size() && (time <= arcPitchYaw.at(arcPitchYaw.size()-1).time))
          logWarning<<"epoch("<<time.dateTimeStr()<<") <= last epoch("<<arcPitchYaw.at(arcPitchYaw.size()-1).time.dateTimeStr()<<")"<<Log::endl;

        StarCameraEpoch epoch;
        epoch.time = time;
        epoch.rotary = rotaryZ(Angle(yaw*1e-6))*rotaryY(Angle(pitch*1e-6));

        arcPitchYaw.push_back(epoch);
      } // for(idEpoch)
    } // for(idFile)

    // =============================================

    logStatus<<"sort epochs"<<Log::endl;
    arcPitchYaw.sort();

    Arc::printStatistics(arcPitchYaw);
    if(arcPitchYaw.size() == 0)
      return;

    if(!fileNameOut.empty())
    {
      logInfo<<"write data to <"<<fileNameOut<<">"<<Log::endl;
      InstrumentFile::write(fileNameOut, arcPitchYaw);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
