/***********************************************/
/**
* @file graceL1b2Accelerometer.cpp
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
This program converts accelerometer data (ACC1B or ACT1B) from the GRACE SDS format into \file{instrument file (ACCELEROMETER)}{instrument}.

Multiple \config{inputfile}s must be given in the correct time order.
The output is one arc of satellite data which can include data gaps.
To split the arc in multiple gap free arcs use \program{InstrumentSynchronize}.

The GRACE SDS format is described in "GRACE Level 1B Data Product User Handbook JPL D-22027"
given at \url{https://podaac-tools.jpl.nasa.gov/drive/files/allData/grace/docs/Handbook_1B_v1.3.pdf}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "fileGrace.h"

/***** CLASS ***********************************/

/** @brief Read GRACE L1B data.
* @ingroup programsConversionGroup */
class GraceL1b2Accelerometer
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GraceL1b2Accelerometer, SINGLEPROCESS, "read GRACE L1B data (ACC1B or ACT1B)", Conversion, Grace, Instrument)

/***********************************************/

void GraceL1b2Accelerometer::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameOutAcc, fileNameOutAng, fileNameOutFlags;
    std::vector<FileName> fileNameIn;

    readConfig(config, "outputfileAccelerometer",        fileNameOutAcc,   Config::OPTIONAL, "", "ACCELEROMETER");
    readConfig(config, "outputfileAngularAccelerometer", fileNameOutAng,   Config::OPTIONAL, "", "ACCELEROMETER");
    readConfig(config, "outputfileFlags",                fileNameOutFlags, Config::OPTIONAL, "", "MISCVALUES(qualflg, acl_res.x, acl_res.y, acl_res.z)");
    readConfig(config, "inputfile",                      fileNameIn,       Config::MUSTSET,  "", "ACC1B or ACT1B");
    if(isCreateSchema(config)) return;

    // =============================================

    logStatus<<"read input files"<<Log::endl;
    Arc arc, arcAngular, arcFlags;
    for(UInt idFile=0; idFile<fileNameIn.size(); idFile++)
    {
      logStatus<<"read file <"<<fileNameIn.at(idFile)<<">"<<Log::endl;
      UInt numberOfRecords;
      FileInGrace file(fileNameIn.at(idFile), numberOfRecords);

      for(UInt idEpoch=0; idEpoch<numberOfRecords; idEpoch++)
      {
        Int32    seconds;
        Char     GRACE_id;
        Vector3d lin_accl, ang_acc, acl_res;
        Byte     qualflg;

        try
        {
          file>>seconds>>GRACE_id>>lin_accl>>ang_acc>>acl_res>>FileInGrace::flag(qualflg);
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

        {
          AccelerometerEpoch epoch;
          epoch.time         = time;
          epoch.acceleration = lin_accl;
          arc.push_back(epoch);
        }

        {
          AccelerometerEpoch epoch;
          epoch.time         = time;
          epoch.acceleration = ang_acc;
          arcAngular.push_back(epoch);
        }

        {
          MiscValuesEpoch epoch(4);
          epoch.time   = time;
          epoch.values = {Double(qualflg), acl_res.x(), acl_res.y(), acl_res.z()};
          arcFlags.push_back(epoch);
        }
      } // for(idEpoch)
    } // for(idFile)

    // =============================================

    logStatus<<"sort epochs"<<Log::endl;
    arc.sort();
    arcAngular.sort();
    arcFlags.sort();

    logStatus<<"eliminate duplicates"<<Log::endl;
    const UInt oldSize = arc.size();
    arc.removeDuplicateEpochs(TRUE/*keepFirst*/);
    arcAngular.removeDuplicateEpochs(TRUE/*keepFirst*/);
    arcFlags.removeDuplicateEpochs(TRUE/*keepFirst*/);
    if(arc.size() < oldSize)
      logInfo<<" "<<oldSize-arc.size()<<" duplicates removed!"<<Log::endl;

    Arc::printStatistics(arc);
    if(arc.size() == 0)
      return;

    if(!fileNameOutAcc.empty())
    {
      logInfo<<"write linear  acceleration to <"<<fileNameOutAcc<<">"<<Log::endl;
      InstrumentFile::write(fileNameOutAcc, arc);
    }
    if(!fileNameOutAng.empty())
    {
      logInfo<<"write angular acceleration to <"<<fileNameOutAng<<">"<<Log::endl;
      InstrumentFile::write(fileNameOutAng, arcAngular);
    }
    if(!fileNameOutFlags.empty())
    {
      logInfo<<"write flags to <"<<fileNameOutFlags<<">"<<Log::endl;
      InstrumentFile::write(fileNameOutFlags, arcFlags);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
