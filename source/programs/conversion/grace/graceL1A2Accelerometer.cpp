/***********************************************/
/**
* @file graceL1A2Accelerometer.cpp
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
This program converts Level-1A accelerometer data to the GROOPS instrument file format.
The GRACE Level-1A format is described in \verb|GRACEiolib.h| given at
\url{http://podaac-tools.jpl.nasa.gov/drive/files/allData/grace/sw/GraceReadSW_L1_2010-03-31.tar.gz}.
Multiple \config{inputfile}s must be given in the correct time order.
The output is one arc of satellite data which can include data gaps.
To split the arc in multiple gap free arcs use \program{InstrumentSynchronize}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "fileGrace.h"

/***** CLASS ***********************************/

/** @brief Read Level-1A GRACE data.
* @ingroup programsConversionGroup */
class GraceL1A2Accelerometer
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(GraceL1A2Accelerometer, SINGLEPROCESS, "read GRACE L1A data", Conversion, Grace, Instrument)

/***********************************************/

void GraceL1A2Accelerometer::run(Config &config)
{
  try
  {
    FileName fileNameAcc, fileNameAccGap, fileNameAngAcc, fileNameHousekeeping;
    std::vector<FileName> fileNameIn;

    readConfig(config, "outputfileAccelerometer",        fileNameAcc,          Config::OPTIONAL, "", "");
    readConfig(config, "outputfileAccelerometerGaps",    fileNameAccGap,       Config::OPTIONAL, "", "");
    readConfig(config, "outputfileAngularAccelerometer", fileNameAngAcc,       Config::OPTIONAL, "", "");
    readConfig(config, "outputfileAccHousekeeping",      fileNameHousekeeping, Config::OPTIONAL, "", "");
    readConfig(config, "inputfile",                      fileNameIn,           Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    // =============================================

    UInt countGaps           = 0;

    logStatus<<"read input files"<<Log::endl;
    Arc arc, arcGap, arcAngAcc, arcHousekeeping;
    for(UInt idFile=0; idFile<fileNameIn.size(); idFile++)
    {
      logStatus<<"read file <"<<fileNameIn.at(idFile)<<">"<<Log::endl;
      UInt numberOfRecords;
      FileInGrace file(fileNameIn.at(idFile), numberOfRecords);

      for(UInt idEpoch=0; idEpoch<numberOfRecords; idEpoch++)
      {
        Int32    seconds, microSeconds, MhzCount;            // seconds, microseconds part, MHz clock count
        Byte     timeRef, GRACE_id, qualityFlag;             // time reference frame (R = Receiver Time, G = GPS time), GRACE satellite ID, data quality flag
        Byte     prodFlag1, prodFlag2, prodFlag3, prodFlag4; // product flag
        Byte     status, TenhzCount;                         // status, 10Hz clock count
        Vector3d acceleration, angularAcceleration;          // linear acceleration, angular acceleration,
        Double   biasVol;                                    // proof mass bias voltage (averaged) (V)
        Float    vd;                                         // amplitude of the AC voltages that operates the position sensors (Vrms)
        Float    x1Out, x2Out, x3Out;                        // displacement of capacitive sensor X1, X2, X3 (m)
        Float    y1Out, y2Out, z1Out;                        // displacement of capacitive sensor Y1, Y2, Y3 (m)
        Float    tesu;                                       // temperature of SU electronics (°C)
        Float    taicu;                                      // temperature of ICU power supply board (°C)
        Float    tisu;                                       // temperature of internal core (°C)
        Float    v15Picu;                                    // ICU reference voltage +15 V
        Float    v15Micu;                                    // ICU reference voltage -15 V
        Float    vr5Picu;                                    // ICU reference voltage + 5 V
        Float    tcicu;                                      // temperature of ICU A/D converter board (°C)
        Float    v15Psu;                                     // SU voltage +15 V
        Float    v15Msu;                                     // SU voltage -15 V
        Float    v48Psu;                                     // SU voltage +48 V
        Float    v48Msu;                                     // SU voltage -48 V
        UInt16   icuBlkNr;                                   // ICU block number

        file>>seconds>>microSeconds;
        file>>timeRef>>GRACE_id>>FileInGrace::flag(qualityFlag)>>FileInGrace::flag(prodFlag1)>>FileInGrace::flag(prodFlag2)>>FileInGrace::flag(prodFlag3)>>FileInGrace::flag(prodFlag4);

        if((Bool(prodFlag4 & (1 << 0)) == 1) && (Bool(prodFlag4 & (1 << 1)) == 1) &&(Bool(prodFlag4 & (1 << 2)) == 1))
          file>>acceleration;
        if((Bool(prodFlag4 & (1 << 3)) == 1) && (Bool(prodFlag4 & (1 << 4)) == 1) &&(Bool(prodFlag4 & (1 << 5)) == 1))
          file>>angularAcceleration;
        if(Bool(prodFlag4 & (1 << 6)) == 1)
          file>>biasVol;
        if(Bool(prodFlag4 & (1 << 7)) == 1)
          file>>vd;

        if(Bool(prodFlag3 & (1 << 0)) == 1)
          file>>x1Out;
        if(Bool(prodFlag3 & (1 << 1)) == 1)
          file>>x2Out;
        if(Bool(prodFlag3 & (1 << 2)) == 1)
          file>>x3Out;
        if(Bool(prodFlag3 & (1 << 3)) == 1)
          file>>y1Out;
        if(Bool(prodFlag3 & (1 << 4)) == 1)
          file>>y2Out;
        if(Bool(prodFlag3 & (1 << 5)) == 1)
          file>>z1Out;
        if(Bool(prodFlag3 & (1 << 6)) == 1)
          file>>tesu;
        if(Bool(prodFlag3 & (1 << 7)) == 1)
          file>>taicu;

        if(Bool(prodFlag2 & (1 << 0)) == 1)
          file>>tisu;
        if(Bool(prodFlag2 & (1 << 1)) == 1)
          file>>v15Picu;
        if(Bool(prodFlag2 & (1 << 2)) == 1)
          file>>v15Micu;
        if(Bool(prodFlag2 & (1 << 3)) == 1)
          file>>vr5Picu;
        if(Bool(prodFlag2 & (1 << 4)) == 1)
          file>>tcicu;
        if(Bool(prodFlag2 & (1 << 5)) == 1)
          file>>v15Psu;
        if(Bool(prodFlag2 & (1 << 6)) == 1)
          file>>v15Msu;
        if(Bool(prodFlag2 & (1 << 7)) == 1)
          file>>v48Psu;

        if(Bool(prodFlag1 & (1 << 0)) == 1)
          file>>v48Msu;
        if(Bool(prodFlag1 & (1 << 1)) == 1)
          file>>status;
        if(Bool(prodFlag1 & (1 << 2)) == 1)
          file>>icuBlkNr;
        if(Bool(prodFlag1 & (1 << 3)) == 1)
          file>>TenhzCount;
        if(Bool(prodFlag1 & (1 << 4)) == 1)
          file>>MhzCount;

        if((Bool(qualityFlag & (1 << 1)) == 0) && (Bool(qualityFlag & (1 << 3)) == 0)) // data with no pulse sync and invalid time tag (?) is removed
        {
          if(microSeconds >= 1'000'000)
          {
            seconds += 1;
            microSeconds -= 1'000'000;
          }

          const Time time = mjd2time(51544.5) + seconds2time(seconds) + seconds2time(microSeconds*1e-6);
          if(arc.size() && (time <= arc.at(arc.size()-1).time))
            logWarning<<"epoch("<<time.dateTimeStr()<<") <= last epoch("<<arc.at(arc.size()-1).time.dateTimeStr()<<")"<<Log::endl;

          // data gaps (> 0.2 seconds)
          if(arc.size() && ((time-arc.at(arc.size()-1).time).seconds()>0.2))
          {
            countGaps++;
            arcGap.push_back(arc.at(arc.size()-1));
          }

          if((Bool(prodFlag4 & (1 << 0)) == 1) && (Bool(prodFlag4 & (1 << 1)) == 1) &&(Bool(prodFlag4 & (1 << 2)) == 1))
          {
            Accelerometer1AEpoch epoch;
            epoch.time         = time;
            epoch.rcvTimeInt   = seconds;
            epoch.rcvTimeFrac  = microSeconds;
            epoch.acceleration = acceleration;
            arc.push_back(epoch);
          }

          if((Bool(prodFlag4 & (1 << 3)) == 1) && (Bool(prodFlag4 & (1 << 4)) == 1) &&(Bool(prodFlag4 & (1 << 5)) == 1))
          {
            Accelerometer1AEpoch epoch;
            epoch.time = time;
            epoch.rcvTimeInt   = seconds;
            epoch.rcvTimeFrac  = microSeconds;
            epoch.acceleration = angularAcceleration;
            arcAngAcc.push_back(epoch);
          }

          if((Bool(prodFlag4 & (1 << 6)) == 1) && (Bool(prodFlag4 & (1 << 7)) == 1) &&
             (Bool(prodFlag3 & (1 << 0)) == 1) && (Bool(prodFlag3 & (1 << 1)) == 1) && (Bool(prodFlag3 & (1 << 2)) == 1) &&
             (Bool(prodFlag3 & (1 << 3)) == 1) && (Bool(prodFlag3 & (1 << 4)) == 1) && (Bool(prodFlag3 & (1 << 5)) == 1) &&
             (Bool(prodFlag3 & (1 << 6)) == 1) && (Bool(prodFlag3 & (1 << 7)) == 1) &&
             (Bool(prodFlag2 & (1 << 0)) == 1) && (Bool(prodFlag2 & (1 << 4)) == 1) && (Bool(prodFlag1 & (1 << 2)) == 1))
          {
            AccHousekeepingEpoch epoch;
            epoch.time = time;
            epoch.biasVoltage  = biasVol;
            epoch.vd           = vd;
            epoch.xOut.x()     = x1Out;
            epoch.xOut.y()     = x2Out;
            epoch.xOut.z()     = x3Out;
            epoch.yOut.x()     = y1Out;
            epoch.yOut.y()     = y2Out;
            epoch.yOut.z()     = z1Out;
            epoch.tempSU       = tesu;
            epoch.tempICU      = taicu;
            epoch.tempCore     = tisu;
            epoch.tempICUConv  = tcicu;
            epoch.blkNrICU     = icuBlkNr;
            arcHousekeeping.push_back(epoch);
          }
        }
      } // for(idEpoch)
    } // for(idFile)

     // =============================================

    logInfo<<"Accelerometer:"<<Log::endl;
    Arc::printStatistics(arc);
    logInfo<<"  gaps (> 0.2 s):  "<<countGaps<<Log::endl;
    if(arc.size() == 0)
      return;

    // remove duplicates
    arc.sort();
    UInt countDuplicates = 0;
    for(UInt idEpoch=1; idEpoch<arc.size(); idEpoch++)
    {
      if(arc.at(idEpoch).time == arc.at(idEpoch-1).time)
      {
        arc.remove(idEpoch--);
        countDuplicates++;
      }
    }
    logInfo<<"  duplicates:      "<<countDuplicates<<Log::endl;

    if(!fileNameAcc.empty())
    {
      logInfo<<"write linear  acceleration to <"<<fileNameAcc<<">"<<Log::endl;
      InstrumentFile::write(fileNameAcc, arc);
    }
    if(!fileNameAccGap.empty())
    {
      logInfo<<"write data gaps to <"<<fileNameAccGap<<">"<<Log::endl;
      InstrumentFile::write(fileNameAccGap, arcGap);
    }

    // =============================================

    logInfo<<"Angular accelerometer:"<<Log::endl;
    Arc::printStatistics(arcAngAcc);

    // remove duplicates
    arc.sort();
    countDuplicates = 0;
    for(UInt idEpoch=1; idEpoch<arcAngAcc.size(); idEpoch++)
    {
      if(arcAngAcc.at(idEpoch).time == arcAngAcc.at(idEpoch-1).time)
      {
        arcAngAcc.remove(idEpoch--);
        countDuplicates++;
      }
    }
    logInfo<<"  duplicates:      "<<countDuplicates<<Log::endl;

    if(!fileNameAngAcc.empty())
    {
      logInfo<<"write angular acceleration to <"<<fileNameAngAcc<<">"<<Log::endl;
      InstrumentFile::write(fileNameAngAcc, arcAngAcc);
    }

    // =============================================

    logInfo<<"Accelerometer housekeeping:"<<Log::endl;
    Arc::printStatistics(arcHousekeeping);
    if(arcHousekeeping.size() == 0)
      return;

    if(!fileNameHousekeeping.empty())
    {
      logInfo<<"write housekeeping data to <"<<fileNameHousekeeping<<">"<<Log::endl;
      InstrumentFile::write(fileNameHousekeeping, arcHousekeeping);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
