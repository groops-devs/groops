/***********************************************/
/**
* @file graceL1b2StarCamera.cpp
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
This program converts orientation data measured by a star camera (SRF to CRF)
from the GRACE SDS format (SCA1B) into \file{instrument file (STARCAMERA)}{instrument}.
For further information see \program{GraceL1b2Accelerometer}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "fileGrace.h"

/***** CLASS ***********************************/

/** @brief Read GRACE L1B data.
* @ingroup programsConversionGroup */
class GraceL1b2StarCamera
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GraceL1b2StarCamera, SINGLEPROCESS, "read GRACE L1B data (SCA1B)", Conversion, Grace, Instrument)

/***********************************************/

void GraceL1b2StarCamera::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName              fileNameSca, FileNameFlags;
    std::vector<FileName> fileNameIn;

    readConfig(config, "outputfileStarCamera",      fileNameSca,   Config::OPTIONAL, "", "");
    readConfig(config, "outputfileStarCameraFlags", FileNameFlags, Config::OPTIONAL, "", "MISCVALUES(sca_id, qual_rss, qualflg)");
    readConfig(config, "inputfile",                 fileNameIn,    Config::MUSTSET,  "", "SCA1B");
    if(isCreateSchema(config)) return;

    // =============================================

    logStatus<<"read input files"<<Log::endl;
    Arc arc, arcFlags;
    for(UInt idFile=0; idFile<fileNameIn.size(); idFile++)
    {
      logStatus<<"read file <"<<fileNameIn.at(idFile)<<">"<<Log::endl;
      UInt numberOfRecords;
      FileInGrace file(fileNameIn.at(idFile), numberOfRecords);

      for(UInt idEpoch=0; idEpoch<numberOfRecords; idEpoch++)
      {
        Int32             seconds;
        Char              GRACE_id;
        FileInGrace::Int8 sca_id;
        Double            quatangle, quaticoeff, quatjcoeff, quatkcoeff, qual_rss;
        Byte              qualflg;

        try
        {
          file>>seconds>>GRACE_id>>sca_id>>quatangle>>quaticoeff>>quatjcoeff>>quatkcoeff>>qual_rss>>FileInGrace::flag(qualflg);
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
          StarCameraEpoch epoch;
          epoch.time   = time;
          epoch.rotary = Rotary3d(Vector{quatangle, quaticoeff, quatjcoeff, quatkcoeff});
          arc.push_back(epoch);
        }

        {
          MiscValuesEpoch epoch(3);
          epoch.time   = time;
          epoch.values = {Double(sca_id), qual_rss, Double(qualflg)};
          arcFlags.push_back(epoch);
        }
      } // for(idEpoch)
    } // for(idFile)

    // =============================================

    logStatus<<"sort epochs"<<Log::endl;
    arc.sort();
    arcFlags.sort();

    logStatus<<"eliminate duplicates"<<Log::endl;
    const UInt oldSize = arc.size();
    arc.removeDuplicateEpochs(TRUE/*keepFirst*/);
    arcFlags.removeDuplicateEpochs(TRUE/*keepFirst*/);
    if(arc.size() < oldSize)
      logInfo<<" "<<oldSize-arc.size()<<" duplicates removed!"<<Log::endl;

    Arc::printStatistics(arc);
    if(arc.size() == 0)
      return;

    if(!fileNameSca.empty())
    {
      logInfo<<"write data to <"<<fileNameSca<<">"<<Log::endl;
      InstrumentFile::write(fileNameSca, arc);
    }

    if(!FileNameFlags.empty())
    {
      logInfo<<"write data to <"<<FileNameFlags<<">"<<Log::endl;
      InstrumentFile::write(FileNameFlags, arcFlags);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
