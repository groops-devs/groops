/***********************************************/
/**
* @file graceL1b2SatelliteTracking.cpp
*
* @brief Read GRACE L1B data.
*
* @author Torsten Mayer-Guerr
* @date 2005-01-19
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program converts low-low satellite data measured by the K-band ranging system
from the GRACE SDS format (KBR1B or LRI1B) into \file{instrument file (SATELLITETRACKING)}{instrument}.
The \config{inputfile}s contain also corrections to antenna offsets
and the so called light time correction. The corrections can be stored in additional files
in the same format as the observations.
If a phase break is found an artificial gap is created.
For further information see \program{GraceL1b2Accelerometer}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "fileGrace.h"

/***** CLASS ***********************************/

/** @brief Read GRACE L1B data.
* @ingroup programsConversionGroup */
class GraceL1b2SatelliteTracking
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GraceL1b2SatelliteTracking, SINGLEPROCESS, "read GRACE L1B data (KBR1B or LRI1B)", Conversion, Grace, Instrument)

/***********************************************/

void GraceL1b2SatelliteTracking::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName              fileNameSst, fileNameAntCentr, fileNameLighttime, fileNameSnr, fileNameIonoCorr;
    std::vector<FileName> fileNameIn;

    readConfig(config, "outputfileSatelliteTracking", fileNameSst,       Config::OPTIONAL, "", "SATELLITETRACKING");
    readConfig(config, "outputfileAntCentr",          fileNameAntCentr,  Config::OPTIONAL, "", "SATELLITETRACKING");
    readConfig(config, "outputfileLighttime",         fileNameLighttime, Config::OPTIONAL, "", "SATELLITETRACKING");
    readConfig(config, "outputfileSNR",               fileNameSnr,       Config::OPTIONAL, "", "MISCVALUES(K_A_SNR, Ka_A_SNR, K_B_SNR, Ka_B_SNR, qualflg)");
    readConfig(config, "outputfileIonoCorr",          fileNameIonoCorr,  Config::OPTIONAL, "", "MISCVALUE");
    readConfig(config, "inputfile",                   fileNameIn,        Config::MUSTSET,  "", "KBR1B or LRI1B");
    if(isCreateSchema(config)) return;

    // =============================================

    logStatus<<"read input files"<<Log::endl;
    SatelliteTrackingArc arc, arcAntCentr, arcLight;
    MiscValuesArc        arcSnr;
    MiscValueArc         arcIonoCorr;
    for(UInt idFile=0; idFile<fileNameIn.size(); idFile++)
    {
      logStatus<<"read file <"<<fileNameIn.at(idFile)<<">"<<Log::endl;
      UInt numberOfRecords;
      FileInGrace file(fileNameIn.at(idFile), numberOfRecords);

      for(UInt idEpoch=0; idEpoch<numberOfRecords; idEpoch++)
      {
        Int32    seconds;
        Double   biased_range, range_rate, range_accl, iono_corr;
        Double   lighttime_corr, lighttime_rate, lighttime_accl;
        Double   ant_centr_corr, ant_centr_rate, ant_centr_accl;
        UInt16   K_A_SNR, Ka_A_SNR, K_B_SNR, Ka_B_SNR;
        Byte     qualflg;

        try
        {
          file>>seconds;
          file>>biased_range>>range_rate>>range_accl>>iono_corr;
          file>>lighttime_corr>>lighttime_rate>>lighttime_accl;
          file>>ant_centr_corr>>ant_centr_rate>>ant_centr_accl;
          file>>K_A_SNR>>Ka_A_SNR>>K_B_SNR>>Ka_B_SNR;
          file>>FileInGrace::flag(qualflg);
        }
        catch(std::exception &/*e*/)
        {
          // GRACE-FO number of records issue
          logWarning<<arc.back().time.dateTimeStr()<<": file ended at "<<idEpoch<<" of "<<numberOfRecords<<" expected records"<<Log::endl;
          break;
        }

        const Time time = mjd2time(51544.5) + seconds2time(seconds);
        if(arc.size() && (time <= arc.back().time))
          logWarning<<"epoch("<<time.dateTimeStr()<<") <= last epoch("<<arc.back().time.dateTimeStr()<<")"<<Log::endl;

        if(((qualflg & 1) == 1) || (arc.size() && ((time-arc.back().time).seconds() <= 5.) && (std::fabs(biased_range-arc.back().range) > 100)))
        {
//           if(arc.size())
//             logWarning<<time.dateTimeStr()<<": skip epoch due to possible phase break, range change = "<<biased_range-arc.back().range<<Log::endl;
//           else
//             logWarning<<time.dateTimeStr()<<": skip epoch due to possible phase break, first epoch in file"<<Log::endl;
          continue;
        }

        {
          SatelliteTrackingEpoch epoch;
          epoch.time              = time;
          epoch.range             = biased_range;
          epoch.rangeRate         = range_rate;
          epoch.rangeAcceleration = range_accl;
          arc.push_back(epoch);
        }

        {
          SatelliteTrackingEpoch epoch;
          epoch.time              = time;
          epoch.range             = lighttime_corr;
          epoch.rangeRate         = lighttime_rate;
          epoch.rangeAcceleration = lighttime_accl;
          arcLight.push_back(epoch);
        }

        {
          SatelliteTrackingEpoch epoch;
          epoch.time              = time;
          epoch.range             = ant_centr_corr;
          epoch.rangeRate         = ant_centr_rate;
          epoch.rangeAcceleration = ant_centr_accl;
          arcAntCentr.push_back(epoch);
        }

        {
          MiscValuesEpoch epoch(5);
          epoch.time = time;
          epoch.values = {Double(K_A_SNR), Double(Ka_A_SNR), Double(K_B_SNR), Double(Ka_B_SNR), Double(qualflg)};
          arcSnr.push_back(epoch);
        }

        {
          MiscValueEpoch epoch;
          epoch.time = time;
          epoch.value = iono_corr;
          arcIonoCorr.push_back(epoch);
        }
      } // for(idEpoch)
    } // for(idFile)

    // =============================================

    logStatus<<"sort epochs"<<Log::endl;
    arc.sort();
    arcAntCentr.sort();
    arcLight.sort();
    arcSnr.sort();
    arcIonoCorr.sort();

    logStatus<<"eliminate duplicates"<<Log::endl;
    const UInt oldSize = arc.size();
    arc.removeDuplicateEpochs(TRUE/*keepFirst*/);
    arcAntCentr.removeDuplicateEpochs(TRUE/*keepFirst*/);
    arcLight.removeDuplicateEpochs(TRUE/*keepFirst*/);
    arcSnr.removeDuplicateEpochs(TRUE/*keepFirst*/);
    if(arc.size() < oldSize)
      logInfo<<" "<<oldSize-arc.size()<<" duplicates removed!"<<Log::endl;

    Arc::printStatistics(arc);
    if(arc.size() == 0)
      return;

    if(!fileNameSst.empty())
    {
      logInfo<<"write data to <"<<fileNameSst<<">"<<Log::endl;
      InstrumentFile::write(fileNameSst, arc);
    }
    if(!fileNameAntCentr.empty())
    {
      logInfo<<"write data to <"<<fileNameAntCentr<<">"<<Log::endl;
      InstrumentFile::write(fileNameAntCentr, arcAntCentr);
    }
    if(!fileNameLighttime.empty())
    {
      logInfo<<"write data to <"<<fileNameLighttime<<">"<<Log::endl;
      InstrumentFile::write(fileNameLighttime, arcLight);
    }
    if(!fileNameSnr.empty())
    {
      logInfo<<"write data to <"<<fileNameSnr<<">"<<Log::endl;
      InstrumentFile::write(fileNameSnr, arcSnr);
    }
    if(!fileNameIonoCorr.empty())
    {
      logInfo<<"write data to <"<<fileNameIonoCorr<<">"<<Log::endl;
      InstrumentFile::write(fileNameIonoCorr, arcIonoCorr);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
