/***********************************************/
/**
* @file graceL1b2GnssReceiver.cpp
*
* @brief Read GRACE L1B data.
*
* @author Torsten Mayer-Guerr
* @date 2012-11-04
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program converts GPS receiver data (phase and pseudo range) data
from the GRACE SDS format (GPS1B or GPS1A) into \file{instrument file (GNSSRECEIVER)}{instrument}.
For further information see \program{GraceL1b2Accelerometer}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "fileGrace.h"

/***** CLASS ***********************************/

/** @brief Read GRACE L1B data.
* @ingroup programsConversionGroup */
class GraceL1b2GnssReceiver
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GraceL1b2GnssReceiver, SINGLEPROCESS, "read GRACE L1B data (GPS1B or GPS1A)", Conversion, Grace, Gnss, Instrument)

/***********************************************/

void GraceL1b2GnssReceiver::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName              fileNameOut;
    std::vector<FileName> fileNameIn;

    readConfig(config, "outputfileGnssReceiver", fileNameOut, Config::MUSTSET,  "", "GNSSRECEIVER");
    readConfig(config, "inputfile",              fileNameIn,  Config::MUSTSET,  "", "GPS1B or GPS1A");
    if(isCreateSchema(config)) return;

    // =============================================

    logStatus<<"read input files"<<Log::endl;
    GnssReceiverArc arc;
    for(UInt idFile=0; idFile<fileNameIn.size(); idFile++)
    {
      logStatus<<"read file <"<<fileNameIn.at(idFile)<<">"<<Log::endl;
      UInt numberOfRecords;
      FileInGrace file(fileNameIn.at(idFile), numberOfRecords);

      for(UInt idEpoch=0; idEpoch<numberOfRecords; idEpoch++)
      {
        Int32             seconds, seconds_frac;
        Char              GRACE_id;
        FileInGrace::Int8 prn_id, ant_id;
        UInt16            prodFlag;
        Byte              qualflg;
        Double            CA_range=NAN_EXPR, L1_range=NAN_EXPR, L2_range=NAN_EXPR;
        Double            CA_phase=NAN_EXPR, L1_phase=NAN_EXPR, L2_phase=NAN_EXPR;
        UInt16            CA_SNR=0, L1_SNR=0, L2_SNR=0, CA_chan=0,L1_chan=0, L2_chan=0;
        Double            K_phase, Ka_phase;                 // K-Band carrier phase, Ka-Band carrier phase
        UInt16            K_SNR, Ka_SNR;                     // K-Band SNR, Ka-Band SNR

        try
        {
          file>>seconds>>seconds_frac>>GRACE_id>>prn_id>>ant_id>>FileInGrace::flag(prodFlag)>>FileInGrace::flag(qualflg);
          if(prodFlag & (1<<0))  file>>CA_range;
          if(prodFlag & (1<<1))  file>>L1_range;
          if(prodFlag & (1<<2))  file>>L2_range;
          if(prodFlag & (1<<3))  file>>CA_phase;
          if(prodFlag & (1<<4))  file>>L1_phase;
          if(prodFlag & (1<<5))  file>>L2_phase;
          if(prodFlag & (1<<6))  file>>CA_SNR;
          if(prodFlag & (1<<7))  file>>L1_SNR;
          if(prodFlag & (1<<8))  file>>L2_SNR;
          if(prodFlag & (1<<9))  file>>CA_chan;
          if(prodFlag & (1<<10)) file>>L1_chan;
          if(prodFlag & (1<<11)) file>>L2_chan;
          if(prodFlag & (1<<12)) file>>K_phase;
          if(prodFlag & (1<<13)) file>>Ka_phase;
          if(prodFlag & (1<<14)) file>>K_SNR;
          if(prodFlag & (1<<15)) file>>Ka_SNR;
        }
        catch(std::exception &/*e*/)
        {
          // GRACE-FO number of records issue
          logWarning<<arc.back().time.dateTimeStr()<<": file ended at "<<idEpoch<<" of "<<numberOfRecords<<" expected records"<<Log::endl;
          break;
        }

        const Time time = mjd2time(51544.5) + seconds2time(seconds) + seconds2time(1e-6*seconds_frac);
        if(arc.size() && (time < arc.back().time))
          throw(Exception("time error: epoch("+time.dateTimeStr()+") < last epoch("+arc.back().time.dateTimeStr()+")"));

        if(ant_id != 0)
        {
          logWarning<<time.dateTimeStr()<<": data from other antennas: "<<(int)ant_id<<Log::endl;
          continue;
        }

        if((prodFlag & 0x3f) != 0x3f)
          logWarning<<time.dateTimeStr()<<" (PRN = "<<prn_id%"%02i), Product: "s<<prodFlag<<Log::endl;

        // start new epoch?
        if((arc.size() == 0) || (arc.back().time != time))
        {
          GnssReceiverEpoch epoch;
          epoch.time    = time;
          epoch.obsType = {GnssType::RANGE + GnssType::L1 + GnssType::C + GnssType::GPS,
                           GnssType::RANGE + GnssType::L1 + GnssType::W + GnssType::GPS,
                           GnssType::RANGE + GnssType::L2 + GnssType::W + GnssType::GPS,
                           GnssType::PHASE + GnssType::L1 + GnssType::C + GnssType::GPS,
                           GnssType::PHASE + GnssType::L1 + GnssType::W + GnssType::GPS,
                           GnssType::PHASE + GnssType::L2 + GnssType::W + GnssType::GPS,
                           GnssType::SNR   + GnssType::L1 + GnssType::C + GnssType::GPS,
                           GnssType::SNR   + GnssType::L1 + GnssType::W + GnssType::GPS,
                           GnssType::SNR   + GnssType::L2 + GnssType::W + GnssType::GPS};
          arc.push_back(epoch);
        }

        arc.back().satellite.push_back(GnssType::GPS + GnssType(prn_id));
        arc.back().observation.push_back(CA_range);
        arc.back().observation.push_back(L1_range);
        arc.back().observation.push_back(L2_range);
        arc.back().observation.push_back(CA_phase);
        arc.back().observation.push_back(L1_phase);
        arc.back().observation.push_back(L2_phase);
        arc.back().observation.push_back(CA_SNR);
        arc.back().observation.push_back(L1_SNR);
        arc.back().observation.push_back(L2_SNR);
      } // for(idEpoch)
    } // for(idFile)

    // =============================================

    Arc::printStatistics(arc);
    if(arc.size() == 0)
      return;

    if(!fileNameOut.empty())
    {
      logInfo<<"write gnss data to <"<<fileNameOut<<">"<<Log::endl;
      InstrumentFile::write(fileNameOut, arc);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
