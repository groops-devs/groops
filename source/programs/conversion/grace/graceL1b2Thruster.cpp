/***********************************************/
/**
* @file graceL1b2Thruster.cpp
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
This program converts thruster data (THR1B or THR1A) from the GRACE SDS format into \file{instrument file (THRUSTER)}{instrument}.
For further information see \program{GraceL1b2Accelerometer}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "fileGrace.h"

/***** CLASS ***********************************/

/** @brief Read GRACE L1B data.
* @ingroup programsConversionGroup */
class GraceL1b2Thruster
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GraceL1b2Thruster, SINGLEPROCESS, "read GRACE L1B data (THR1B or THR1A)", Conversion, Grace, Instrument)

/***********************************************/

void GraceL1b2Thruster::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName              fileNameOut;
    std::vector<FileName> fileNameIn;

    readConfig(config, "outputfileThruster", fileNameOut, Config::MUSTSET,  "", "THRUSTER");
    readConfig(config, "inputfile",          fileNameIn,  Config::MUSTSET,  "", "THR1B or THR1A");
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
        constexpr UInt MAXTHRSTRS = 14;
        Int32    seconds, time_frac;
        Char     time_ref, GRACE_id;
        UInt32   thrust_count[MAXTHRSTRS], on_time[MAXTHRSTRS], accum_dur[MAXTHRSTRS];
        Byte     qualflg;

        try
        {
          file>>seconds>>time_frac>>time_ref>>GRACE_id;
          for(UInt k=0; k<MAXTHRSTRS; k++)
            file>>thrust_count[k];
          for(UInt k=0; k<MAXTHRSTRS; k++)
            file>>on_time[k];
          for(UInt k=0; k<MAXTHRSTRS; k++)
            file>>accum_dur[k];
          file>>FileInGrace::flag(qualflg);
        }
        catch(std::exception &/*e*/)
        {
          // GRACE-FO number of records issue
          logWarning<<arc.back().time.dateTimeStr()<<": file ended at "<<idEpoch<<" of "<<numberOfRecords<<" expected records"<<Log::endl;
          break;
        }

        const Time time = mjd2time(51544.5) + seconds2time(seconds) + seconds2time(1e-6*time_frac);
        if(arc.size() && (time <= arc.back().time))
          logWarning<<"epoch("<<time.dateTimeStr()<<") <= last epoch("<<arc.back().time.dateTimeStr()<<")"<<Log::endl;

//         if(time_ref != 'G')
//         {
//           logWarning<<time.dateTimeStr()<<" time is not GPS time: '"<<time_ref<<"' (ignored)"<<Log::endl;
//           continue;
//         }

        ThrusterEpoch epoch;
        epoch.time     = time;
        epoch.onTime1  = on_time[0];
        epoch.onTime2  = on_time[1];
        epoch.onTime3  = on_time[2];
        epoch.onTime4  = on_time[3];
        epoch.onTime5  = on_time[4];
        epoch.onTime6  = on_time[5];
        epoch.onTime7  = on_time[6];
        epoch.onTime8  = on_time[7];
        epoch.onTime9  = on_time[8];
        epoch.onTime10 = on_time[9];
        epoch.onTime11 = on_time[10];
        epoch.onTime12 = on_time[11];
        epoch.onTime13 = on_time[12];
        epoch.onTime14 = on_time[13];
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
