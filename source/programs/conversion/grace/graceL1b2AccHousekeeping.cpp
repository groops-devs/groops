/***********************************************/
/**
* @file graceL1b2AccHousekeeping.cpp
*
* @brief Read GRACE L1B data.
*
* @author Beate Klinger
* @date 2015-10-08
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program converts ACC housekeeping data (AHK1B or AHK1A) from the GRACE SDS format into \file{instrument file (ACCHOUSEKEEPING)}{instrument}.
For further information see \program{GraceL1b2Accelerometer}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "fileGrace.h"

/***** CLASS ***********************************/

/** @brief Read GRACE L1B data.
* @ingroup programsConversionGroup */
class GraceL1b2AccHousekeeping
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GraceL1b2AccHousekeeping, SINGLEPROCESS, "read GRACE L1B data (AHK1B or AHK1A)", Conversion, Grace, Instrument)

/***********************************************/

void GraceL1b2AccHousekeeping::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName              fileNameOut;
    std::vector<FileName> fileNameIn;

    readConfig(config, "outputfileAccHousekeeping", fileNameOut, Config::MUSTSET,  "", "ACCHOUSEKEEPING");
    readConfig(config, "inputfile",                 fileNameIn,  Config::MUSTSET,  "", "AHK1B or AHK1A");
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
        Int32             seconds, time_frac;
        Char              time_ref, GRACE_id;
        Byte              qualflg;
        UInt32            prodFlag;
        Vector3d          acceleration, angularAcceleration; // linear acceleration, angular acceleration,
        Double            biasVol=NAN_EXPR;                  // proof mass bias voltage (averaged) (V)
        Float             vd=NAN_EXPR;                       // amplitude of the AC voltages that operates the position sensors (Vrms)
        Float             x1Out=NAN_EXPR, x2Out=NAN_EXPR, x3Out=NAN_EXPR;   // displacement of capacitive sensor X1, X2, X3 (m)
        Float             y1Out=NAN_EXPR, y2Out=NAN_EXPR, z1Out=NAN_EXPR;   // displacement of capacitive sensor Y1, Y2, Y3 (m)
        Float             tesu =NAN_EXPR;                    // temperature of SU electronics (°C)
        Float             taicu=NAN_EXPR;                    // temperature of ICU power supply board (°C)
        Float             tisu =NAN_EXPR;                    // temperature of internal core (°C)
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
          file>>seconds>>time_frac>>time_ref>>GRACE_id>>FileInGrace::flag(qualflg)>>FileInGrace::flag(prodFlag);
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

//         if(!((prodFlag & (1 << 14)) && (prodFlag & (1 << 15)) && (prodFlag & (1 << 16)) && (prodFlag & (1 << 20))))
//           continue;

        AccHousekeepingEpoch epoch;
        epoch.time        = mjd2time(51544.5) + seconds2time(seconds) + seconds2time(1e-6*time_frac);
        epoch.biasVoltage = biasVol;
        epoch.vd          = vd;
        epoch.xOut        = Vector3d(x1Out, x2Out, x3Out);
        epoch.yOut        = Vector3d(y1Out, y2Out, z1Out);
        epoch.tempSU      = tesu;
        epoch.tempICU     = taicu;
        epoch.tempCore    = tisu;
        epoch.tempICUConv = tcicu;
        epoch.blkNrICU    = icuBlkNr;
        arc.push_back(epoch);
      } // for(idEpoch)
    } // for(idFile)

    // =============================================

    logStatus<<"sort epochs"<<Log::endl;
    arc.sort();

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
