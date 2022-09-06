/***********************************************/
/**
* @file instrumentAccelerometerApplyEstimatedParameters.cpp
*
* @brief Apply estimated acceleration parameters to an accelerometer file.
*
* @author Torsten Mayer-Guerr
* @date 2016-02-16
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program evaluates estimated satellite parameters and writes the result to an accelerometer file.
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

/** @brief Apply estimated acceleration parameters to an accelerometer file
* @ingroup programsGroup */
class InstrumentAccelerometerApplyEstimatedParameters
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(InstrumentAccelerometerApplyEstimatedParameters, PARALLEL, "apply estimated acceleration parameters to an accelerometer file.", Instrument)
GROOPS_RENAMED_PROGRAM(InstrumentAccelerometerEstimatedParameter, InstrumentAccelerometerApplyEstimatedParameters, date2time(2022, 8, 8))

/***********************************************/

void InstrumentAccelerometerApplyEstimatedParameters::run(Config &config, Parallel::CommunicatorPtr comm)
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
    Double                         factor;

    renameDeprecatedConfig(config, "satelliteModel", "inputfileSatelliteModel",     date2time(2020, 8, 19));
    renameDeprecatedConfig(config, "parameter",      "parametrizationAcceleration", date2time(2020, 6,  3));

    readConfig(config, "outputfileAccelerometer",     fileNameOutAccelerometer, Config::MUSTSET,  "", "");
    readConfig(config, "inputfileSatelliteModel",     fileNameSatellite,        Config::OPTIONAL, "{groopsDataDir}/satelliteModel/", "satellite macro model");
    readConfig(config, "inputfileOrbit",              fileNameOrbit,            Config::OPTIONAL, "",  "");
    readConfig(config, "inputfileStarCamera",         fileNameStarCamera,       Config::OPTIONAL, "",  "");
    readConfig(config, "inputfileAccelerometer",      fileNameAccelerometer,    Config::OPTIONAL, "",  "add non-gravitational forces in satellite reference frame");
    readConfig(config, "earthRotation",               earthRotation,            Config::OPTIONAL, "",  "");
    readConfig(config, "ephemerides",                 ephemerides,              Config::OPTIONAL, "jpl", "may be needed by parametrizationAcceleration");
    readConfig(config, "parametrizationAcceleration", parameterAcceleration,    Config::MUSTSET,  "",  "orbit force parameters");
    readConfig(config, "inputfileParameter",          fileNameSolution,         Config::MUSTSET,  "",  "estimated orbit force parameters");
    readConfig(config, "indexStart",                  indexStart,               Config::DEFAULT,  "0", "position in the solution vector");
    readConfig(config, "rightSide",                   rightSide,                Config::DEFAULT,  "0", "if solution contains several right hand sides, select one");
    readConfig(config, "factor",                      factor,                   Config::DEFAULT, "1.0", "the result is multiplied by this factor");
    if(isCreateSchema(config)) return;

    // open and test instrument files
    // ------------------------------
    InstrumentFile orbitFile(fileNameOrbit);
    InstrumentFile starCameraFile(fileNameStarCamera);
    InstrumentFile accelerometerFile(fileNameAccelerometer);
    InstrumentFile::checkArcCount({orbitFile, starCameraFile, accelerometerFile});
    UInt arcCount = orbitFile.arcCount();
    if(!arcCount) arcCount = starCameraFile.arcCount();
    if(!arcCount) arcCount = accelerometerFile.arcCount();

    SatelliteModelPtr satellite;
    if(!fileNameSatellite.empty())
      readFileSatelliteModel(fileNameSatellite, satellite);

    Matrix mx;
    readFileMatrix(fileNameSolution, mx);
    Vector parameter = mx.slice((indexStart>=0) ? UInt(indexStart) : mx.rows()-UInt(-indexStart),
                                (rightSide>=0)  ? UInt(rightSide) : mx.columns()-UInt(-rightSide),
                                parameterAcceleration->parameterCount() + arcCount*parameterAcceleration->parameterCountArc(),
                                1);

    // Accelerometer-Daten erzeugen
    // -----------------------
    logStatus<<"computing accelerations"<<Log::endl;
    std::vector<Arc> arcList(arcCount);
    Parallel::forEach(arcList, [&](UInt arcNo)
    {
      OrbitArc         orbit         = orbitFile.readArc(arcNo);
      StarCameraArc    starCamera    = starCameraFile.readArc(arcNo);
      AccelerometerArc accelerometer = accelerometerFile.readArc(arcNo);
      Arc::checkSynchronized({orbit, starCamera, accelerometer});
      std::vector<Time> times = orbit.times();
      if(!times.size()) times = starCamera.times();
      if(!times.size()) times = accelerometer.times();

      parameterAcceleration->setIntervalArc(times.front(), times.back()+medianSampling(times));
      const UInt countA = parameterAcceleration->parameterCount();
      const UInt countB = parameterAcceleration->parameterCountArc();

      AccelerometerArc accArc;
      for(UInt i=0; i<times.size(); i++)
      {
        Rotary3d rotEarth, rotSat;
        Vector3d position, velocity;
        if(earthRotation)     rotEarth = earthRotation->rotaryMatrix(times.at(i));
        if(starCamera.size()) rotSat   = starCamera.at(i).rotary;
        if(orbit.size())      position = orbit.at(i).position;
        if(orbit.size())      velocity = orbit.at(i).velocity;

        Matrix A(3, countA);
        Matrix B(3, countB);
        parameterAcceleration->compute(satellite, times.at(i), position, velocity, rotSat, rotEarth, ephemerides, A, B);
        Vector g(3);
        if(countA) matMult(factor, A, parameter.row(0, countA), g);
        if(countB) matMult(factor, B, parameter.row(arcNo*countB+countA, countB), g);

        AccelerometerEpoch epoch;
        epoch.time         = times.at(i);
        epoch.acceleration = rotSat.inverseRotate(rotEarth.inverseRotate(Vector3d(g)));
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
