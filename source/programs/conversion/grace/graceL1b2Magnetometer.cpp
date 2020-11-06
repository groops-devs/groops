/***********************************************/
/**
* @file graceL1b2Magnetometer.cpp
*
* @brief Read GRACE L1B data.
*
* @author Beate Klinger
* @date 2014-05-26
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program converts magnetometer data from the GRACE SDS format into \file{instrument file (MAGNETOMETER)}{instrument}.
For further information see \program{GraceL1b2Accelerometer}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "fileGrace.h"

/***** CLASS ***********************************/

/** @brief Read GRACE L1B data.
* @ingroup programsConversionGroup */
class GraceL1b2Magnetometer
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(GraceL1b2Magnetometer, SINGLEPROCESS, "read GRACE L1B data", Conversion, Grace, Instrument)

/***********************************************/

void GraceL1b2Magnetometer::run(Config &config)
{
  try
  {
    FileName              fileNameOut;
    std::vector<FileName> fileNameIn;

    readConfig(config, "outputfileMagnetometer", fileNameOut, Config::MUSTSET,  "", "");
    readConfig(config, "inputfile",              fileNameIn,  Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    // =============================================

    logStatus<<"read input files"<<Log::endl;
    Arc arc;
    for(UInt idFile=0; idFile<fileNameIn.size(); idFile++)
    {
      logStatus<<"read file <"<<fileNameIn.at(idFile)<<">"<<Log::endl;
      UInt numberOfRecords;
      FileInGrace file(fileNameIn.at(idFile), numberOfRecords);

      for(UInt idEpoch=0; idEpoch<numberOfRecords; idEpoch++)
      {
        Int32    seconds, time_frac;
        Byte     time_ref, GRACE_id;
        Float    MfvX_RAW, MfvY_RAW, MfvZ_RAW;
        Float    torque1A, torque2A, torque3A, torque1B, torque2B, torque3B;
        Float    MF_BCalX, MF_BCalY, MF_BCalZ, torque_cal;
        Byte     qualflg;

        file>>seconds>>time_frac>>time_ref>>GRACE_id;
        file>>MfvX_RAW>>MfvY_RAW>>MfvZ_RAW;
        file>>torque1A>>torque2A>>torque3A>>torque1B>>torque2B>>torque3B;
        file>>MF_BCalX>>MF_BCalY>>MF_BCalZ>>torque_cal>>FileInGrace::flag(qualflg);

        const Time time = mjd2time(51544.5) + seconds2time(seconds) + seconds2time(1e-6*time_frac);
        if(arc.size() && (time <= arc.at(arc.size()-1).time))
          logWarning<<"epoch("<<time.dateTimeStr()<<") <= last epoch("<<arc.at(arc.size()-1).time.dateTimeStr()<<")"<<Log::endl;

        if(time_ref != 'G')
        {
          logWarning<<time.dateTimeStr()<<" time is not GPS time"<<Log::endl;
          continue;
        }

        MagnetometerEpoch epoch;
        epoch.time                     = time;
        epoch.magneticField            = Vector3d(MfvX_RAW, MfvY_RAW, MfvZ_RAW);
        epoch.torquerA                 = Vector3d(torque1A, torque2A, torque3A);
        epoch.torquerB                 = Vector3d(torque1B, torque2B, torque3B);
        epoch.magneticFieldCalibration = Vector3d(MF_BCalX, MF_BCalY, MF_BCalZ);
        epoch.torquerCalibration       = torque_cal;
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
