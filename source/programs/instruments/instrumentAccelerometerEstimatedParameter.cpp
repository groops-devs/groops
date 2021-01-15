/***********************************************/
/**
* @file instrumentAccelerometerEstimatedParameter.cpp
*
* @brief Estimated satellite parameter as accelerometer data.
*
* @author Torsten Mayer-Guerr
* @date 2016-02-16
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program evaluates estimated satellite parameters and write the result as accelerometer file.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"
#include "files/fileSatelliteModel.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/ephemerides/ephemerides.h"
#include "classes/parametrizationAcceleration/parametrizationAcceleration.h"

/***** CLASS ***********************************/

/** @brief Estimated satellite parameter as accelerometer data.
* @ingroup programsGroup */
class InstrumentAccelerometerEstimatedParameter
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(InstrumentAccelerometerEstimatedParameter, PARALLEL, "Estimated satellite parameter as accelerometer data", Instrument)

/***********************************************/

void InstrumentAccelerometerEstimatedParameter::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName                       fileNameOutAccelerometer;
    FileName                       fileNameSatellite;
    FileName                       fileNameOrbit, fileNameStarCamera, fileNameAccelerometer;
    EarthRotationPtr               earthRotation;
    EphemeridesPtr                 ephemerides;
    ParametrizationAccelerationPtr parameterAcceleration;
    FileName                       fileNameSolution;
    Int                            rightSide, indexStart;

    renameDeprecatedConfig(config, "satelliteModel", "inputfileSatelliteModel",     date2time(2020, 8, 19));
    renameDeprecatedConfig(config, "parameter",      "parametrizationAcceleration", date2time(2020, 6,  3));

    readConfig(config, "outputfileAccelerometer",     fileNameOutAccelerometer, Config::MUSTSET,  "", "");
    readConfig(config, "inputfileSatelliteModel",     fileNameSatellite,        Config::OPTIONAL, "{groopsDataDir}/satelliteModel/", "satellite macro model");
    readConfig(config, "inputfileOrbit",              fileNameOrbit,            Config::MUSTSET,  "",  "");
    readConfig(config, "inputfileStarCamera",         fileNameStarCamera,       Config::MUSTSET,  "",  "");
    readConfig(config, "inputfileAccelerometer",      fileNameAccelerometer,    Config::OPTIONAL, "",  "add non-gravitational forces in satellite reference frame");
    readConfig(config, "earthRotation",               earthRotation,            Config::MUSTSET,  "",  "");
    readConfig(config, "ephemerides",                 ephemerides,              Config::OPTIONAL, "jpl", "may be needed by parametrizationAcceleration");
    readConfig(config, "parametrizationAcceleration", parameterAcceleration,    Config::MUSTSET,  "",  "orbit force parameters");
    readConfig(config, "inputfileParameter",          fileNameSolution,         Config::MUSTSET,  "",  "estimated orbit force parameters");
    readConfig(config, "indexStart",                  indexStart,               Config::DEFAULT,  "0", "position in the solution vector");
    readConfig(config, "rightSide",                   rightSide,                Config::DEFAULT,  "0", "if solution contains several right hand sides, select one");
    if(isCreateSchema(config)) return;

    // open and test instrument files
    // ------------------------------
    InstrumentFile orbitFile(fileNameOrbit);
    InstrumentFile starCameraFile(fileNameStarCamera);
    InstrumentFile accelerometerFile(fileNameAccelerometer);
    InstrumentFile::checkArcCount({orbitFile, starCameraFile, accelerometerFile});

    SatelliteModelPtr satellite;
    if(!fileNameSatellite.empty())
      readFileSatelliteModel(fileNameSatellite, satellite);

    Matrix mx;
    readFileMatrix(fileNameSolution, mx);
    Vector parameter = mx.slice((indexStart>=0) ? UInt(indexStart) : mx.rows()-UInt(-indexStart),
                                (rightSide>=0)  ? UInt(rightSide) : mx.columns()-UInt(-rightSide),
                                parameterAcceleration->parameterCount() + orbitFile.arcCount()*parameterAcceleration->parameterCountArc(),
                                1);

    // Accelerometer-Daten erzeugen
    // -----------------------
    logStatus<<"computing accelerations"<<Log::endl;
    std::vector<Arc> arcList(orbitFile.arcCount());
    Parallel::forEach(arcList, [&](UInt arcNo)
    {
      OrbitArc         orbit         = orbitFile.readArc(arcNo);
      StarCameraArc    starCamera    = starCameraFile.readArc(arcNo);
      AccelerometerArc accelerometer = accelerometerFile.readArc(arcNo);
      Arc::checkSynchronized({orbit, starCamera, accelerometer});
      const UInt epochCount = orbit.size();

      parameterAcceleration->setIntervalArc(orbit.at(0).time, orbit.back().time+medianSampling(orbit.times()));
      const UInt countA = parameterAcceleration->parameterCount();
      const UInt countB = parameterAcceleration->parameterCountArc();

      AccelerometerArc accArc;
      for(UInt i=0; i<epochCount; i++)
      {
        const Rotary3d rotEarth = earthRotation->rotaryMatrix(orbit.at(i).time);
        Matrix A(3, countA);
        Matrix B(3, countB);
        parameterAcceleration->compute(satellite, orbit.at(i).time, orbit.at(i).position, orbit.at(i).velocity,
                                       starCamera.at(i).rotary, rotEarth, ephemerides, A, B);
        Vector g2(3);
        if(countA) matMult(1., A, parameter.row(0, countA), g2);
        if(countB) matMult(1., B, parameter.row(arcNo*countB+countA, countB), g2);

        AccelerometerEpoch epoch;
        epoch.time         = orbit.at(i).time;
        epoch.acceleration = starCamera.at(i).rotary.inverseRotate(rotEarth.inverseRotate(Vector3d(g2(0), g2(1), g2(2))));
        if(accelerometer.size())
          epoch.acceleration += accelerometer.at(i).acceleration;
        accArc.push_back(epoch);
      }

      return accArc;
    }, comm);

    // write result
    // ------------
    if(Parallel::isMaster(comm))
    {
      logStatus<<"write accelerometer data to file <"<<fileNameOutAccelerometer<<">"<<Log::endl;
      InstrumentFile::write(fileNameOutAccelerometer, arcList);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
