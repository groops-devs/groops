/***********************************************/
/**
* @file graceL1A2Temperature.cpp
*
* @brief Read GRACE L1A data.
*
* @author Beate Klinger
* @date 2017-03-09
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program converts Level-1A temperature measurments to the GROOPS instrument file format.
The GRACE Level-1A format is described in GRACE given at \url{http://podaac-tools.jpl.nasa.gov/drive/files/allData/grace/sw/GraceReadSW_L1_2010-03-31.tar.gz}.
Multiple \config{inputfile}s must be given in the correct time order.
The output is one arc of satellite data which can include data gaps.
To split the arc in multiple gap free arcs use \program{InstrumentSynchronize}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "fileGrace.h"

/***** CLASS ***********************************/

/** @brief Read GRACE L1A data.
* @ingroup programsConversionGroup */
class GraceL1A2Temperature
{
  void readFile(const FileName &fileName, MiscValuesArc &arcTemp);

public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GraceL1A2Temperature, SINGLEPROCESS, "read GRACE L1A data", Conversion, Grace, Instrument)

/***********************************************/

void GraceL1A2Temperature::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameTemp;
    std::vector<FileName> fileNameIn;

    readConfig(config, "outputfileTemperature",   fileNameTemp,   Config::OPTIONAL, "", "");
    readConfig(config, "inputfile",               fileNameIn,        Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    // =============================================

    logStatus<<"read input files"<<Log::endl;
    Arc arc;
    for(UInt idFile=0; idFile<fileNameIn.size(); idFile++)
    {
      logStatus<<"read file '"<<fileNameIn.at(idFile)<<"'"<<Log::endl;
      UInt numberOfRecords;
      FileInGrace file(fileNameIn.at(idFile), numberOfRecords);

      for(UInt idEpoch=0; idEpoch<numberOfRecords; idEpoch++)
      {
        Int32 seconds, microSeconds;           // seconds, microseconds part
        Byte  timeRef, GRACE_id, qualityFlag;  // time reference frame (R = Receiver Time, G = GPS time), GRACE satellite ID, data quality flag
        Float TempMEPNegY;                     // I/F support structure to MEP -y  [°C]
        Float TempMEPPosY;                     // I/F support structure to MEP +y  [°C]
        Float TempMEPm;                        // I/F support structure to MEP mid [°C]
        Float TempICU;                         // ICU Temperature Reference Point  [°C]
        Float TempICURed;                      // ICU Tempeartuer Reference Point (redundant) [°C]
        Float TempACCNegZ;                     // ACC thermal cage -z [°C]
        Float TempACCPosZ;                     // ACC thermal cage +z [°C]
        Float TempCFRPPosX, TempCFRPPosXRed;   // CFRP Frame at +x I/F to Baseplate (+redundant) [°C]
        Float TempCFRPNegX, TempCFRPNegXRed;   // CFRP Frame at -x I/F to Baseplate (+redundant) [°C]
        Float TempCFRPNegY, TempCFRPNegYRed;   // CFRP Frame at -y I/F to Baseplate (+redundant) [°C]
        Float TempACCSen;                      // Harness to ACC sensor [°C]
        Float TempICUSpec;                     // ICU special [°C]
        Float TempMWANegY,  TempMWANegYOff;    // MWA -y Temp. Ref. Point, 2 out 3 (+nominally off) [°C]
        Float TempMWAPosY,  TempMWAPosYOff;    // MWA +y Temp. Ref. Point, 2 out 3 (+nominally off) [°C]
        Float TempHornPosX, TempHornPosXRed;   // Horn aperture +x (+redundant) [°C]
        Float TempHornPl,   TempHornPlRed;     // Horn / platform I/F (+redundant) [°C]
        Float TempHWMANegY;                    // Harnass to MWA electronics -y [°C]
        Float TempHWMAPosY;                    // Harnass to MWA electronics +y [°C]
        Float TempRFSamp;                      // RF-Sampling unit [°C]
        Float TempUSONegY,  TempUSONegYRed;    // USO Temp. Ref. Point -y (+redundant) [°C]
        Float TempUSOPosY,  TempUSOPosYRed;    // USO Temp. Ref. Point +y (+redundant) [°C]

        file>>seconds>>microSeconds>>timeRef>>GRACE_id>>TempMEPNegY>>TempMEPPosY>>TempMEPm>>TempICU>>TempICURed>>TempACCNegZ>>TempACCPosZ;
        file>>TempCFRPPosX>>TempCFRPPosXRed>>TempCFRPNegX>>TempCFRPNegXRed>>TempCFRPNegY>>TempCFRPNegYRed>>TempACCSen>>TempICUSpec;
        file>>TempMWANegY>>TempMWANegYOff>>TempMWAPosY>>TempMWAPosYOff>>TempHornPosX>>TempHornPosXRed>>TempHornPl>>TempHornPlRed;
        file>>TempHWMANegY>>TempHWMAPosY>>TempRFSamp>>TempUSONegY>>TempUSONegYRed>>TempUSOPosY>>TempUSOPosYRed>>FileInGrace::flag(qualityFlag);

        {
          MiscValuesEpoch epoch(21);
          epoch.time = mjd2time(51544.5) + seconds2time(seconds) + seconds2time(microSeconds*1e-6);
          epoch.values(0)  = TempACCNegZ;
          epoch.values(1)  = TempACCPosZ;
          epoch.values(2)  = TempACCSen;
          epoch.values(3)  = TempCFRPPosX;
          epoch.values(4)  = TempCFRPNegX;
          epoch.values(5)  = TempCFRPNegY;
          epoch.values(6)  = TempCFRPNegYRed;
          epoch.values(7)  = TempMEPNegY;
          epoch.values(8)  = TempMEPPosY;
          epoch.values(9)  = TempMEPm;
          epoch.values(10) = TempICUSpec;
          epoch.values(11) = TempICU;
          epoch.values(12) = TempMWANegY;
          epoch.values(13) = TempMWAPosY;
          epoch.values(14) = TempHWMANegY;
          epoch.values(15) = TempHWMAPosY;
          epoch.values(16) = TempRFSamp;
          epoch.values(17) = TempUSONegY;
          epoch.values(18) = TempUSOPosY;
          epoch.values(19) = TempHornPosX;
          epoch.values(20) = TempHornPl;
          arc.push_back(epoch);
        }
      }// for(idEpoch)
    } // for(idFile)

    // =============================================

    Arc::printStatistics(arc);
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

    if(!fileNameTemp.empty())
    {
      logInfo<<"write temperature data to <"<<fileNameTemp<<">"<<Log::endl;
      InstrumentFile::write(fileNameTemp, arc);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
