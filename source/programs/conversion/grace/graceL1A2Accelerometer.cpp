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
This program converts Level-1A accelerometer data (ACC1A) to the GROOPS instrument file format.
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
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GraceL1A2Accelerometer, SINGLEPROCESS, "read GRACE L1A data (ACC1A)", Conversion, Grace, Instrument)

/***********************************************/

void GraceL1A2Accelerometer::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameOutAcc, fileNameOutAng;
    std::vector<FileName> fileNameIn;

    readConfig(config, "outputfileAccelerometer",        fileNameOutAcc, Config::OPTIONAL, "", "ACCELEROMETER1A");
    readConfig(config, "outputfileAngularAccelerometer", fileNameOutAng, Config::OPTIONAL, "", "ACCELEROMETER1A");
    readConfig(config, "inputfile",                      fileNameIn,     Config::MUSTSET,  "", "ACC1A");
    if(isCreateSchema(config)) return;

    // =============================================

    logStatus<<"read input files"<<Log::endl;
    Arc arc, arcAngAcc;
    for(UInt idFile=0; idFile<fileNameIn.size(); idFile++)
    {
      logStatus<<"read file <"<<fileNameIn.at(idFile)<<">"<<Log::endl;
      UInt numberOfRecords;
      FileInGrace file(fileNameIn.at(idFile), numberOfRecords);

      for(UInt idEpoch=0; idEpoch<numberOfRecords; idEpoch++)
      {
        Int32             seconds, microSeconds;
        Char              time_ref, GRACE_id;
        Byte              qualflg;
        UInt32            prodFlag;
        Vector3d          acceleration, angularAcceleration; // linear acceleration, angular acceleration,
        Double            biasVol=NAN_EXPR;                  // proof mass bias voltage (averaged) (V)
        Float             vd=NAN_EXPR;                       // amplitude of the AC voltages that operates the position sensors (Vrms)
        Float             x1Out, x2Out, x3Out;               // displacement of capacitive sensor X1, X2, X3 (m)
        Float             y1Out, y2Out, z1Out;               // displacement of capacitive sensor Y1, Y2, Y3 (m)
        Float             tesu;                              // temperature of SU electronics (°C)
        Float             taicu;                             // temperature of ICU power supply board (°C)
        Float             tisu;                              // temperature of internal core (°C)
        Float             v15Picu;                           // ICU reference voltage +15 V
        Float             v15Micu;                           // ICU reference voltage -15 V
        Float             vr5Picu;                           // ICU reference voltage + 5 V
        Float             tcicu;                             // temperature of ICU A/D converter board (°C)
        Float             v15Psu;                            // SU voltage +15 V
        Float             v15Msu;                            // SU voltage -15 V
        Float             v48Psu;                            // SU voltage +48 V
        Float             v48Msu;                            // SU voltage -48 V
        Byte              status;                            // status
        UInt16            icuBlkNr=0;                        // ICU block number
        FileInGrace::Int8 PPS_source;                        // 10Hz clock count
        Int32             sync_quality_index;
        UInt32            status_flag;

        try
        {
          file>>seconds>>microSeconds>>time_ref>>GRACE_id>>FileInGrace::flag(qualflg)>>FileInGrace::flag(prodFlag);
          if(prodFlag & (1 <<  0))  file>>acceleration.x();
          if(prodFlag & (1 <<  1))  file>>acceleration.y();
          if(prodFlag & (1 <<  2))  file>>acceleration.z();
          if(prodFlag & (1 <<  3))  file>>angularAcceleration.x();
          if(prodFlag & (1 <<  4))  file>>angularAcceleration.y();
          if(prodFlag & (1 <<  5))  file>>angularAcceleration.z();
          if(prodFlag & (1 <<  6))  file>>biasVol;
          if(prodFlag & (1 <<  7))  file>>vd;
          if(prodFlag & (1 <<  8))  file>>x1Out;
          if(prodFlag & (1 <<  9))  file>>x2Out;
          if(prodFlag & (1 << 10))  file>>x3Out;
          if(prodFlag & (1 << 11))  file>>y1Out;
          if(prodFlag & (1 << 12))  file>>y2Out;
          if(prodFlag & (1 << 13))  file>>z1Out;
          if(prodFlag & (1 << 14))  file>>tesu;
          if(prodFlag & (1 << 15))  file>>taicu;
          if(prodFlag & (1 << 16))  file>>tisu;
          if(prodFlag & (1 << 17))  file>>v15Picu;
          if(prodFlag & (1 << 18))  file>>v15Micu;
          if(prodFlag & (1 << 19))  file>>vr5Picu;
          if(prodFlag & (1 << 20))  file>>tcicu;
          if(prodFlag & (1 << 21))  file>>v15Psu;
          if(prodFlag & (1 << 22))  file>>v15Msu;
          if(prodFlag & (1 << 23))  file>>v48Psu;
          if(prodFlag & (1 << 24))  file>>v48Msu;
          if(prodFlag & (1 << 25))  file>>FileInGrace::flag(status);
          if(prodFlag & (1 << 26))  file>>icuBlkNr;
          if(prodFlag & (1 << 27))  file>>PPS_source;
          if(prodFlag & (1 << 28))  file>>sync_quality_index;
          if(prodFlag & (1 << 29))  file>>FileInGrace::flag(status_flag);
        }
        catch(std::exception &/*e*/)
        {
          // GRACE-FO number of records issue
          logWarning<<arc.back().time.dateTimeStr()<<": file ended at "<<idEpoch<<" of "<<numberOfRecords<<" expected records"<<Log::endl;
          break;
        }

        if((qualflg & (1 << 1)) || (qualflg & (1 << 3))) // data with no pulse sync and invalid time tag (?) are removed
          continue;

        while(microSeconds >= 1'000'000)
        {
          seconds += 1;
          microSeconds -= 1'000'000;
        }

        const Time time = mjd2time(51544.5) + seconds2time(seconds) + seconds2time(microSeconds*1e-6);
        if(arc.size() && (time <= arc.back().time))
          logWarning<<"epoch("<<time.dateTimeStr()<<") <= last epoch("<<arc.back().time.dateTimeStr()<<")"<<Log::endl;

        if((prodFlag & (1<<0)) && (prodFlag & (1<<1)) && (prodFlag & (1<<2)))
        {
          Accelerometer1AEpoch epoch;
          epoch.time         = time;
          epoch.rcvTimeInt   = seconds;
          epoch.rcvTimeFrac  = microSeconds;
          epoch.acceleration = acceleration;
          arc.push_back(epoch);
        }

        if((prodFlag & (1<<3)) && (prodFlag & (1<<4)) && (prodFlag & (1<<5)))
        {
          Accelerometer1AEpoch epoch;
          epoch.time = time;
          epoch.rcvTimeInt   = seconds;
          epoch.rcvTimeFrac  = microSeconds;
          epoch.acceleration = angularAcceleration;
          arcAngAcc.push_back(epoch);
        }
      } // for(idEpoch)
    } // for(idFile)

     // =============================================

    logStatus<<"sort epochs"<<Log::endl;
    arc.sort();
    arcAngAcc.sort();

    logStatus<<"eliminate duplicates"<<Log::endl;
    const UInt oldSize = arc.size();
    arc.removeDuplicateEpochs(TRUE/*keepFirst*/);
    arcAngAcc.removeDuplicateEpochs(TRUE/*keepFirst*/);
    if(arc.size() < oldSize)
      logInfo<<" "<<oldSize-arc.size()<<" duplicates removed!"<<Log::endl;

    Arc::printStatistics(arc);
    if(arc.size() == 0)
      return;

    if(!fileNameOutAcc.empty())
    {
      logInfo<<"write linear  acceleration to <"<<fileNameOutAcc<<">"<<Log::endl;
      InstrumentFile::write(fileNameOutAcc, arc);
    }
    if(!fileNameOutAng.empty())
    {
      logInfo<<"write angular acceleration to <"<<fileNameOutAng<<">"<<Log::endl;
      InstrumentFile::write(fileNameOutAng, arcAngAcc);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
