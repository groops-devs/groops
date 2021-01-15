/***********************************************/
/**
* @file graceL1b2Time.cpp
*
* @brief Read GRACE L1B data.
*
* @author Beate Klinger
* @date 2017-02-02
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program converts time data from the GRACE SDS format into \file{instrument file (MISCVALUES)}{instrument}.
For further information see \program{GraceL1b2Accelerometer}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "fileGrace.h"

/***** CLASS ***********************************/

/** @brief Read GRACE L1B data.
* @ingroup programsConversionGroup */
class GraceL1b2Time
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GraceL1b2Time, SINGLEPROCESS, "read GRACE L1B data", Conversion, Grace, Instrument)

/***********************************************/

void GraceL1b2Time::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName              fileNameOut;
    std::vector<FileName> fileNameIn;

    readConfig(config, "outputfileTime", fileNameOut, Config::MUSTSET,  "", "");
    readConfig(config, "inputfile",      fileNameIn,  Config::MUSTSET,  "", "");
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
        Int32 obdh_time;
        Byte  GRACE_id;
        Int32 TS_suppid, gpstime_intg, gpstime_frac, first_icu_blknr, final_icu_blknr;
        Byte  qualflg;

        file>>obdh_time>>GRACE_id>>TS_suppid>>gpstime_intg>>gpstime_frac>>first_icu_blknr>>final_icu_blknr>>FileInGrace::flag(qualflg);

        const Time time = mjd2time(51544.5) + seconds2time(gpstime_intg);// + seconds2time(1e-6*gpstime_frac);

        if((first_icu_blknr == -1.0) || (final_icu_blknr == -1.0)) // no ACC measurements
          continue;

        MiscValuesEpoch epoch(7);
        epoch.time   = time;
        epoch.values = {Double(obdh_time), Double(TS_suppid), Double(gpstime_intg), Double(gpstime_frac), Double(first_icu_blknr), Double(final_icu_blknr), Double(qualflg)};
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
