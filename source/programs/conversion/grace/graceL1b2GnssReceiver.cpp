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
from the GRACE SDS format into \file{instrument file (GNSSRECEIVER)}{instrument}.
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
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(GraceL1b2GnssReceiver, SINGLEPROCESS, "read GRACE L1B data", Conversion, Grace, Gnss, Instrument)

/***********************************************/

void GraceL1b2GnssReceiver::run(Config &config)
{
  try
  {
    FileName              fileNameOut;
    std::vector<FileName> fileNameIn;

    readConfig(config, "outputfileGnssReceiver", fileNameOut, Config::MUSTSET,  "", "");
    readConfig(config, "inputfile",              fileNameIn,  Config::MUSTSET,  "", "");
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
        Int32    seconds, seconds_frac;
        Byte     GRACE_id, prn_id, ant_id, qualflg;
        UInt16   prod_flag;
        Double   CA_range, L1_range, L2_range, CA_phase, L1_phase, L2_phase;
        UInt16   CA_SNR,   L1_SNR,   L2_SNR,   CA_chan,  L1_chan,  L2_chan;

        file>>seconds>>seconds_frac>>GRACE_id>>prn_id>>ant_id>>FileInGrace::flag(prod_flag)>>FileInGrace::flag(qualflg);
        file>>CA_range>>L1_range>>L2_range>>CA_phase>>L1_phase>>L2_phase;
        file>>CA_SNR>>L1_SNR>>L2_SNR>>CA_chan>>L1_chan>>L2_chan;

        const Time time = mjd2time(51544.5) + seconds2time(seconds) + seconds2time(1e-6*seconds_frac);
        if(arc.size() && (time < arc.at(arc.size()-1).time))
          throw(Exception("time error: epoch("+time.dateTimeStr()+") < last epoch("+arc.at(arc.size()-1).time.dateTimeStr()+")"));

        if(ant_id != 0)
        {
          logWarning<<time.dateTimeStr()<<": data from other antennas: "<<(int)ant_id<<Log::endl;
          continue;
        }

        if((prod_flag & 0x3f) != 0x3f)
          logWarning<<time.dateTimeStr()<<" (PRN = "<<prn_id%"%02i), Product: "s<<prod_flag<<Log::endl;

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

        arc.at(arc.size()-1).satellite.push_back(GnssType::GPS + GnssType(prn_id));
        arc.at(arc.size()-1).observation.push_back(CA_range);
        arc.at(arc.size()-1).observation.push_back(L1_range);
        arc.at(arc.size()-1).observation.push_back(L2_range);
        arc.at(arc.size()-1).observation.push_back(CA_phase);
        arc.at(arc.size()-1).observation.push_back(L1_phase);
        arc.at(arc.size()-1).observation.push_back(L2_phase);
        arc.at(arc.size()-1).observation.push_back(CA_SNR);
        arc.at(arc.size()-1).observation.push_back(L1_SNR);
        arc.at(arc.size()-1).observation.push_back(L2_SNR);
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
