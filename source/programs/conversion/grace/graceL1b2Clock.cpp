/***********************************************/
/**
* @file graceL1b2Clock.cpp
*
* @brief Read GRACE L1B data.
*
* @author Beate Klinger
* @date 2017-01-04
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program converts clock data from the GRACE SDS format into \file{instrument file (CLOCK)}{instrument}.
For further information see \program{GraceL1b2Accelerometer}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "fileGrace.h"

/***** CLASS ***********************************/

/** @brief Read GRACE L1B data.
* @ingroup programsConversionGroup */
class GraceL1b2Clock
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(GraceL1b2Clock, SINGLEPROCESS, "read GRACE L1B data", Conversion, Grace, Instrument)

/***********************************************/

void GraceL1b2Clock::run(Config &config)
{
  try
  {
    FileName              fileNameOut;
    std::vector<FileName> fileNameIn;

    readConfig(config, "outputfileClock", fileNameOut, Config::MUSTSET,  "", "");
    readConfig(config, "inputfile",       fileNameIn,  Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    // =============================================

    logStatus<<"read input files"<<Log::endl;
    Arc arc;
    for(UInt idFile=0; idFile<fileNameIn.size(); idFile++)
    {
      logStatus<<"read file <"<<fileNameIn.at(idFile)<<">"<<Log::endl;
      UInt numberOfRecords;
      FileInGrace file(fileNameIn.at(idFile), numberOfRecords);

      Bool dummy = FALSE;
      for(UInt idEpoch=0; idEpoch<numberOfRecords; idEpoch++)
      {
        Int32    seconds;
        Byte     GRACE_id, clock_id;
        Double   eps_time, eps_err, eps_drift, drift_err;
        Byte     qualflg;

        file>>seconds>>GRACE_id>>clock_id>>eps_time>>eps_err>>eps_drift>>drift_err>>FileInGrace::flag(qualflg);

        const Time time = mjd2time(51544.5) + seconds2time(seconds);
        if(arc.size() && (time <= arc.at(arc.size()-1).time))
          logWarning<<"epoch("<<time.dateTimeStr()<<") <= last epoch("<<arc.at(arc.size()-1).time.dateTimeStr()<<")"<<Log::endl;

        UInt quality = 0;
        // invalid time periods
        if(Bool(qualflg & (1 << 0)) == 1)
          quality = 1;
        if(Bool(qualflg & (1 << 1)) == 1)
          quality = 2;

        // boundaries - days
        if((idEpoch==0) && (quality == 2))
          quality = 0;
        if((arc.size()) && (idEpoch==1) && (time == arc.at(arc.size()-1).time))
          continue;

        if((arc.size()) && (idEpoch==numberOfRecords-1) && (quality == 1) && (time == arc.at(arc.size()-1).time))
          continue;

        ClockEpoch epoch;
        epoch.time        = time;
        epoch.rcvTime     = seconds;
        epoch.epsTime     = eps_time;
        epoch.epsError    = eps_err;
        epoch.epsDrift    = eps_drift;
        epoch.driftError  = drift_err;
        epoch.qualityFlag = quality;
        //arc.push_back(epoch);
        //epoch.values = {Double(clock_id), eps_time, eps_err, eps_drift, drift_err, Double(qualflg)};

         // delete invalid periods
        if(quality == 1)
        {
          dummy = TRUE;
          arc.push_back(epoch);
          continue;
        }
        if((dummy) && (quality != 2))
          continue;
        if((dummy) && (quality == 2))
        {
          dummy = FALSE;
          if((arc.size()) && (epoch.time == arc.at(arc.size()-1).time) && (arc.at(arc.size()-1).data().at(1) != 0))
            continue;
          if((arc.size()) && (epoch.time == arc.at(arc.size()-1).time) && (arc.at(arc.size()-1).data().at(1) == 0))
          {
            arc.remove(arc.size()-1,1);
            if(epoch.epsTime != 0)
            {
              arc.push_back(epoch);
              continue;
            }
          }
        }
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
