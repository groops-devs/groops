/***********************************************/
/**
* @file graceL1A2SatelliteTracking.cpp
*
* @brief Read GRACE L1A data.
*
* @author Beate Klinger
* @date 2018-01-29
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program converts Level-1A satellite tracking data to the GROOPS instrument file format.
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
class GraceL1A2SatelliteTracking
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(GraceL1A2SatelliteTracking, SINGLEPROCESS, "read GRACE L1A data", Conversion, Grace, Instrument)

/***********************************************/

void GraceL1A2SatelliteTracking::run(Config &config)
{
  try
  {
    FileName fileNameSst;
    std::vector<FileName> fileNameIn;

    readConfig(config, "outputfileSatelliteTracking", fileNameSst,    Config::OPTIONAL, "", "");
    readConfig(config, "inputfile", fileNameIn, Config::MUSTSET,  "", "");
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
        Int32    seconds, microSeconds;             // receiver time, microseconds part
        Byte     GRACE_id, prn_id, ant_id;          // GRACE satellite ID, GPS PRN id or GRACE id, KBR antenna id
        Byte     prodFlag1, prodFlag2, qualityFlag; // product flag, quality flag
        Double   CA_range, L1_range, L2_range;
        Double   CA_phase, L1_phase, L2_phase;
        UInt16   CA_SNR, L1_SNR, L2_SNR;
        UInt16   CA_chan, L1_chan, L2_chan;
        Double   K_phase, Ka_phase;                 // K-Band carrier phase, Ka-Band carrier phase
        UInt16   K_SNR, Ka_SNR;                     // K-Band SNR, Ka-Band SNR

        file>>seconds>>microSeconds>>GRACE_id>>prn_id>>ant_id>>FileInGrace::flag(prodFlag1)>>FileInGrace::flag(prodFlag2)>>FileInGrace::flag(qualityFlag);
        if(prodFlag2 & (1<<0)) file>>CA_range;
        if(prodFlag2 & (1<<1)) file>>L1_range;
        if(prodFlag2 & (1<<2)) file>>L2_range;
        if(prodFlag2 & (1<<3)) file>>CA_phase;
        if(prodFlag2 & (1<<4)) file>>L1_phase;
        if(prodFlag2 & (1<<5)) file>>L2_phase;
        if(prodFlag2 & (1<<6)) file>>CA_SNR;
        if(prodFlag2 & (1<<7)) file>>L1_SNR;
        if(prodFlag1 & (1<<0)) file>>L2_SNR;
        if(prodFlag1 & (1<<1)) file>>CA_chan;
        if(prodFlag1 & (1<<2)) file>>L1_chan;
        if(prodFlag1 & (1<<3)) file>>L2_chan;
        if(prodFlag1 & (1<<4)) file>>K_phase;
        if(prodFlag1 & (1<<5)) file>>Ka_phase;
        if(prodFlag1 & (1<<6)) file>>K_SNR;
        if(prodFlag1 & (1<<7)) file>>Ka_SNR;

        //logInfo<<(mjd2time(51544.5) + seconds2time(seconds) + seconds2time(microSeconds*1e-6)).dateTimeStr()<<"  "<<GRACE_id<<"  "<<(int)prn_id<<"  "<<(int)ant_id<<"  "<<K_phase<<"  "<<Ka_phase<<"  "<<K_SNR<<"  "<<Ka_SNR<<Log::endl;

        if(microSeconds >= 1'000'000)
          {
            seconds += 1;
            microSeconds -= 1'000'000;
          }

        if((ant_id == 11) || (ant_id == -11))
        {
          MiscValuesEpoch epoch(7);
          epoch.time      = mjd2time(51544.5) + seconds2time(seconds) + seconds2time(microSeconds*1e-6);
          epoch.values(0) = seconds;
          epoch.values(1) = microSeconds;
          epoch.values(2) = ant_id;
          epoch.values(3) = K_phase;
          epoch.values(4) = Ka_phase;
          epoch.values(5) = K_SNR;
          epoch.values(6) = Ka_SNR;
          arc.push_back(epoch);
        }

        /*if((Bool(qualityFlag & (1 << 1)) == 0) && (Bool(qualityFlag & (1 << 3)) == 0)) // data with no pulse sync and invalid time tag (?) is removed
        {
          if(microSeconds >= 1e6)
          {
            seconds += 1;
            microSeconds -= 1e6;
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

        }*/
      } // for(idEpoch)
    } // for(idFile)

    // =============================================

    arc.sort();
    Arc::printStatistics(arc);
    if(arc.size() == 0)
      return;

    if(!fileNameSst.empty())
    {
      logInfo<<"write data to <"<<fileNameSst<<">"<<Log::endl;
      InstrumentFile::write(fileNameSst, arc);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
