/***********************************************/
/**
* @file instrumentStarCamera2RollPitchYaw.cpp
*
* @brief Compute roll, pitch, yaw angles.
*
* @author Torsten Mayer-Guerr
* @date 2017-06-14
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Compute roll, pitch, yaw angles from \configFile{inputfileStarCamera}{instrument} data.
Optional the angles are computed relative
to a \configFile{inputfileStarCameraReference}{instrument}.

See also \program{SimulateStarCamera}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Compute roll, pitch, yaw angles.
* @ingroup programsGroup */
class InstrumentStarCamera2RollPitchYaw
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(InstrumentStarCamera2RollPitchYaw, PARALLEL, "Compute roll, pitch, yaw angles", Instrument)

/***********************************************/

void InstrumentStarCamera2RollPitchYaw::run(Config &config)
{
  try
  {
    FileName fileNameOut;
    FileName fileNameStarCamera, fileNameReference;

    readConfig(config, "outputfileInstrument",         fileNameOut,        Config::MUSTSET,  "", "roll, pitch, yaw [rad], VECTOR3D");
    readConfig(config, "inputfileStarCamera",          fileNameStarCamera, Config::MUSTSET,  "", "");
    readConfig(config, "inputfileStarCameraReference", fileNameReference,  Config::OPTIONAL, "", "nominal orientation");
    if(isCreateSchema(config)) return;

    logStatus<<"read star camera data and compute roll, pitch, yaw"<<Log::endl;
    InstrumentFile  starCamera1File(fileNameStarCamera);
    InstrumentFile  starCamera2File(fileNameReference);
    InstrumentFile::checkArcCount({starCamera1File, starCamera2File});

    std::vector<Arc> arcList(starCamera1File.arcCount());
    Parallel::forEach(arcList, [&](UInt arcNo)
    {
      StarCameraArc starCamera1 = starCamera1File.readArc(arcNo);
      StarCameraArc starCamera2 = starCamera2File.readArc(arcNo);
      Arc::checkSynchronized({starCamera1, starCamera2});

      Vector3dArc arc;
      for(UInt i=0; i<starCamera1.size(); i++)
      {
        if(starCamera2.size())
          starCamera1.at(i).rotary = inverse(starCamera2.at(i).rotary) * starCamera1.at(i).rotary;

        Angle roll, pitch, yaw;
        inverse(starCamera1.at(i).rotary).cardan(roll, pitch, yaw);

        Vector3dEpoch epoch;
        epoch.time     = starCamera1.at(i).time;
        epoch.vector3d = Vector3d(roll, pitch, yaw);
        arc.push_back(epoch);
      }
      return arc;
    }); // forEach

    if(Parallel::isMaster())
    {
      logStatus<<"write angles to file <"<<fileNameOut<<">"<<Log::endl;
      InstrumentFile::write(fileNameOut, arcList);
      Arc::printStatistics(arcList);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
