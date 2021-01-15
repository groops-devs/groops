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
This program converts ACC housekeeping data from the GRACE SDS format into \file{instrument file (ACCHOUSEKEEPING)}{instrument}.
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

GROOPS_REGISTER_PROGRAM(GraceL1b2AccHousekeeping, SINGLEPROCESS, "read GRACE L1B data", Conversion, Grace, Instrument)

/***********************************************/

void GraceL1b2AccHousekeeping::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName              fileNameOut;
    std::vector<FileName> fileNameIn;

    readConfig(config, "outputfileAccHousekeeping", fileNameOut, Config::MUSTSET,  "", "");
    readConfig(config, "inputfile",                 fileNameIn,  Config::MUSTSET,  "", "");
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
        Byte     qualflg, prodFlag1, prodFlag2, prodFlag3, prodFlag4; // product flag
        Byte     status;                // status
        Double   biasVol;               // proof mass bias voltage (averaged) (V)
        Float    vd;                    // amplitude of the AC voltages that operates the position sensors (Vrms)
        Float    x1Out, x2Out, x3Out;   // displacement of capacitive sensor X1, X2, X3 (m)
        Float    y1Out, y2Out, z1Out;   // displacement of capacitive sensor Y1, Y2, Y3 (m)
        Float    tesu;                  // temperature of SU electronics (°C)
        Float    taicu;                 // temperature of ICU power supply board (°C)
        Float    tisu;                  // temperature of internal core (°C)
        Float    v15Picu;               // ICU reference voltage +15 V
        Float    v15Micu;               // ICU reference voltage -15 V
        Float    vr5Picu;               // ICU reference voltage + 5 V
        Float    tcicu;                 // temperature of ICU A/D converter board (°C)
        Float    v15Psu;                // SU voltage +15 V
        Float    v15Msu;                // SU voltage -15 V
        Float    v48Psu;                // SU voltage +48 V
        Float    v48Msu;                // SU voltage -48 V
        UInt16   icuBlkNr;              // ICU block number

        file>>seconds>>time_frac>>time_ref>>GRACE_id>>FileInGrace::flag(qualflg);
        file>>FileInGrace::flag(prodFlag1)>>FileInGrace::flag(prodFlag2)>>FileInGrace::flag(prodFlag3)>>FileInGrace::flag(prodFlag4);
        if(prodFlag4 & (1<<6)) file>>biasVol;
        if(prodFlag4 & (1<<7)) file>>vd;
        if(prodFlag3 & (1<<0)) file>>x1Out;
        if(prodFlag3 & (1<<1)) file>>x2Out;
        if(prodFlag3 & (1<<2)) file>>x3Out;
        if(prodFlag3 & (1<<3)) file>>y1Out;
        if(prodFlag3 & (1<<4)) file>>y2Out;
        if(prodFlag3 & (1<<5)) file>>z1Out;
        if(prodFlag3 & (1<<6)) file>>tesu;
        if(prodFlag3 & (1<<7)) file>>taicu;
        if(prodFlag2 & (1<<0)) file>>tisu;
        if(prodFlag2 & (1<<1)) file>>v15Picu;
        if(prodFlag2 & (1<<2)) file>>v15Micu;
        if(prodFlag2 & (1<<3)) file>>vr5Picu;
        if(prodFlag2 & (1<<4)) file>>tcicu;
        if(prodFlag2 & (1<<5)) file>>v15Psu;
        if(prodFlag2 & (1<<6)) file>>v15Msu;
        if(prodFlag2 & (1<<7)) file>>v48Psu;
        if(prodFlag1 & (1<<0)) file>>v48Msu;
        if(prodFlag1 & (1<<1)) file>>FileInGrace::flag(status);
        if(prodFlag1 & (1<<2)) file>>icuBlkNr;

        const Time time = mjd2time(51544.5) + seconds2time(seconds) + seconds2time(1e-6*time_frac);

        if(time_ref != 'G')
        {
          logWarning<<time.dateTimeStr()<<" time is not GPS time"<<Log::endl;
          continue;
        }

        if(!((prodFlag3&(1<<6)) && (prodFlag3&(1<<7)) && (prodFlag2&(1<<0)) && (prodFlag2&(1<<4))))
          continue;

        AccHousekeepingEpoch epoch;
        epoch.time        = time;
        epoch.biasVoltage = biasVol;
        epoch.vd          = vd;
        epoch.xOut.x()    = x1Out;
        epoch.xOut.y()    = x2Out;
        epoch.xOut.z()    = x3Out;
        epoch.yOut.x()    = y1Out;
        epoch.yOut.y()    = y2Out;
        epoch.yOut.z()    = z1Out;
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
