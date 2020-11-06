/***********************************************/
/**
* @file noiseStarCamera.cpp
*
* @brief Add noise to star camera observations.
*
* @author Torsten Mayer-Guerr
* @author Matthias Ellmer
* @date 2011-05-24
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program adds noise to rotation observations. The noise is computed via a pseudo random sequence.
See \configClass{noiseGenerator}{noiseGeneratorType} for details on noise options.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "classes/noiseGenerator/noiseGenerator.h"

/***** CLASS ***********************************/

/** @brief Add white noise to quaternion observations.
  * @ingroup programsGroup */
class NoiseStarCamera
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(NoiseStarCamera, PARALLEL, "add noise to rotation observations", Simulation, Noise, Instrument)

/***********************************************/

void NoiseStarCamera::run(Config &config)
{
  try
  {
    FileName inName, outName,  outNameCovariance;
    NoiseGeneratorPtr noiseRoll, noisePitch, noiseYaw;

    readConfig(config, "outputfileStarCamera", outName,    Config::MUSTSET, "", "");
    readConfig(config, "inputfileStarCamera",  inName,     Config::MUSTSET, "", "");
    readConfig(config, "noiseRoll",            noiseRoll,  Config::MUSTSET, "", "[rad]");
    readConfig(config, "noisePitch",           noisePitch, Config::MUSTSET, "", "[rad]");
    readConfig(config, "noiseYaw",             noiseYaw,   Config::MUSTSET, "", "[rad]");
    if(isCreateSchema(config)) return;

    logStatus<<"add noise to star Camera data"<<Log::endl;
    InstrumentFile scaFile(inName);
    UInt arcCount = scaFile.arcCount();
    std::vector<Arc> arcList(arcCount);

    Parallel::forEach(arcList, [&](UInt arcNo)
    {
      StarCameraArc scaArc = scaFile.readArc(arcNo);

      const UInt   posCount = scaArc.size();
      const Vector eRoll    = noiseRoll->noise(posCount);
      const Vector ePitch   = noisePitch->noise(posCount);
      const Vector eYaw     = noiseYaw->noise(posCount);

      StarCameraArc arc;
      for(UInt i=0; i<posCount; i++)
      {
        Angle roll, pitch, yaw;
        scaArc.at(i).rotary.cardan(roll, pitch, yaw);

        Angle nRoll  = Angle(roll  + eRoll(i));
        Angle nPitch = Angle(pitch + ePitch(i));
        Angle nYaw   = Angle(yaw   + eYaw(i));

        StarCameraEpoch epoch;
        epoch.time   = scaArc.at(i).time;
        epoch.rotary = rotaryZ(nYaw) * rotaryY(nPitch) * rotaryX(nRoll);
        arc.push_back(epoch);
      }
      return arc;
    }); // forEach

    logStatus<<"write star camera data to file <"<<outName<<">"<<Log::endl;
    if(Parallel::isMaster())
      InstrumentFile::write(outName, arcList);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
