/***********************************************/
/**
* @file graceL1b2Orbit.cpp
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
This program converts the reduced dynamical orbit
from the GRACE/GRACE-FO SDS format (GNV1B, GNI1B) into \file{instrument file (ORBIT)}{instrument}.

When GNV1B is used, the orbit can be rotated from the terrestrial reference frame (TRF) transformed into the celestial reference frame (CRF) by
specifying \configClass{earthRotation}{earthRotationType}.

For further information see \program{GraceL1b2Accelerometer}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "classes/earthRotation/earthRotation.h"
#include "fileGrace.h"

/***** CLASS ***********************************/

/** @brief Read GRACE/GRACE-FO L1B data.
* @ingroup programsConversionGroup */
class GraceL1b2Orbit
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GraceL1b2Orbit, SINGLEPROCESS, "read GRACE/GRACE-FO L1B data (GNV1B, GNI1B)", Conversion, Grace, Orbit, Instrument)

/***********************************************/

void GraceL1b2Orbit::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName              fileNameOut;
    std::vector<FileName> fileNameIn;
    EarthRotationPtr      earthRotation;

    readConfig(config, "outputfileOrbit", fileNameOut,   Config::MUSTSET,  "", "");
    readConfig(config, "earthRotation",   earthRotation, Config::OPTIONAL, "file", "to rotate GNV1B into CRF");
    readConfig(config, "inputfile",       fileNameIn,    Config::MUSTSET,  "", "GNV1B/GNI1B");
    if(isCreateSchema(config)) return;

    // =============================================

    logStatus<<"read input files"<<Log::endl;
    OrbitArc arc;
    for(UInt idFile=0; idFile<fileNameIn.size(); idFile++)
    {
      logStatus<<"read file <"<<fileNameIn.at(idFile)<<">"<<Log::endl;
      UInt numberOfRecords;
      FileInGrace file(fileNameIn.at(idFile), numberOfRecords);

      for(UInt idEpoch=0; idEpoch<numberOfRecords; idEpoch++)
      {
        Int32    seconds;
        Char     GRACE_id, coord_ref;
        Vector3d pos, pos_err, vel, vel_err;
        Byte     qualflg;

        try
        {
          file>>seconds>>GRACE_id>>coord_ref>>pos>>pos_err>>vel>>vel_err>>FileInGrace::flag(qualflg);
        }
        catch(std::exception &/*e*/)
        {
          // GRACE-FO number of records issue
          logWarning<<arc.back().time.dateTimeStr()<<": file ended at "<<idEpoch<<" of "<<numberOfRecords<<" expected records"<<Log::endl;
          break;
        }

        const Time time = mjd2time(51544.5) + seconds2time(seconds);
        if(arc.size() && (time <= arc.back().time))
          continue;

        if(coord_ref == 'E' && !earthRotation)
          throw(Exception(time.dateTimeStr()+": orbit not in CRF. Must specify earthRotation."));

        OrbitEpoch epoch;
        epoch.time     = time;
        epoch.position = pos;
        epoch.velocity = vel;
        arc.push_back(epoch);
      } // for(idEpoch)
    } // for(idFile)

    // =============================================

    if(earthRotation)
    {
      logStatus<<"rotation from TRF to CRF"<<Log::endl;
      Single::forEach(arc.size(), [&](UInt i)
      {
        const Rotary3d rotation = inverse(earthRotation->rotaryMatrix(arc.at(i).time));
        arc.at(i).position = rotation.rotate(arc.at(i).position);
        if(arc.at(i).velocity.r() > 0)
          arc.at(i).velocity = rotation.rotate(arc.at(i).velocity) + crossProduct(earthRotation->rotaryAxis(arc.at(i).time), arc.at(i).position);
      });
    }

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
