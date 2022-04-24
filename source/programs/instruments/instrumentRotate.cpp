/***********************************************/
/**
* @file instrumentRotate.cpp
*
* @brief Rotate instrument data into a new reference frame.
*
* @author Torsten Mayer-Guerr
* @date 2017-06-15
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program rotates \file{instrument data}{instrument} into a new reference frame
(using \configFile{inputfileStarCamera}{instrument}).
The rotation is usually done from satellite frame into inertial frame.
To apply Earth rotation use \program{InstrumentEarthRotation}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Rotate instrument data into a new reference frame.
* @ingroup programsGroup */
class InstrumentRotate
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(InstrumentRotate, PARALLEL, "Rotate instrument data into a new reference frame", Instrument)

/***********************************************/

void InstrumentRotate::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName fileNameOut, fileNameIn, fileNameStarCamera;
    Bool     inverseRotate;

    readConfig(config, "outputfileInstrument",fileNameOut,        Config::MUSTSET,  "",  "");
    readConfig(config, "inputfileInstrument", fileNameIn,         Config::MUSTSET,  "",  "");
    readConfig(config, "inputfileStarCamera", fileNameStarCamera, Config::MUSTSET,  "",  "");
    readConfig(config, "inverseRotate",       inverseRotate,      Config::DEFAULT,  "1", "");
    if(isCreateSchema(config)) return;

    // read orbit and rotate
    // ---------------------
    logStatus<<"read and rotate data <"<<fileNameIn<<">"<<Log::endl;
    InstrumentFile instrumentFile(fileNameIn);
    InstrumentFile starCameraFile(fileNameStarCamera);
    InstrumentFile::checkArcCount({instrumentFile, starCameraFile});

    std::vector<Arc> arcList(instrumentFile.arcCount());
    Parallel::forEach(arcList, [&](UInt arcNo)
    {
      Arc           arc        = instrumentFile.readArc(arcNo);
      StarCameraArc starCamera = starCameraFile.readArc(arcNo);
      Arc::checkSynchronized({arc, starCamera});

      for(UInt i=0; i<arc.size(); i++)
      {
        Rotary3d rot = starCamera.at(i).rotary;
        if(inverseRotate)
          rot = inverse(rot);

        if(arc.getType() == Epoch::ORBIT)
        {
          auto &epoch = dynamic_cast<OrbitEpoch&>(arc.at(i));
          epoch.position     = rot.rotate(epoch.position);
          epoch.velocity     = rot.rotate(epoch.velocity);
          epoch.acceleration = rot.rotate(epoch.acceleration);
        }
        else if(arc.getType() == Epoch::ACCELEROMETER)
        {
          auto &epoch = dynamic_cast<AccelerometerEpoch&>(arc.at(i));
          epoch.acceleration = rot.rotate(epoch.acceleration);
        }
        else if(arc.getType() == Epoch::GRADIOMETER)
        {
          auto &epoch = dynamic_cast<GradiometerEpoch&>(arc.at(i));
          epoch.gravityGradient = rot.rotate(epoch.gravityGradient);
        }
        else if(arc.getType() == Epoch::COVARIANCE3D)
        {
          auto &epoch = dynamic_cast<Covariance3dEpoch&>(arc.at(i));
          epoch.covariance = rot.rotate(epoch.covariance);
        }
        else if(arc.getType() == Epoch::VECTOR3D)
        {
          auto &epoch = dynamic_cast<Vector3dEpoch&>(arc.at(i));
          epoch.vector3d = rot.rotate(epoch.vector3d);
        }
        else
          throw(Exception("rotation for "+arc.getTypeName()+" not implemented"));
      }
      return arc;
    }, comm); // forEach

    if(Parallel::isMaster(comm))
    {
      logStatus<<"write data to file <"<<fileNameOut<<">"<<Log::endl;
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
