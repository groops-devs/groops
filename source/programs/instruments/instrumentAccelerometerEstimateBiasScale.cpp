/***********************************************/
/**
* @file instrumentAccelerometerEstimateBiasScale.cpp
*
* @brief Estimate accelerometer bias and scale.
*
* @author Beate Klinger
* @date 2014-11-27
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program calibrates \configFile{inputfileAccelerometer}{instrument} with respect to
simulated accelerometer data, see \program{SimulateAccelerometer}.
The parameters \configFile{outputfileSolution}{matrix}
of \configClass{parametrizationAcceleration}{parametrizationAccelerationType}
are estimated and the effect is reduced to calibrate the \file{accelerometer data}{instrument}.

If \configFile{inputfileThruster}{instrument} is given, the corresponding epochs
(within \config{marginThruster}) are not used for the parameter estimation,
but the accelerometer epochs are still calibrated afterwards.
An arbitrary instrument file is allowed here.

The \configFile{inputfileOrbit}{instrument}, \configFile{inputfileStarCamera}{instrument},
\configClass{earthRotation}{earthRotationType}, \configClass{ephemerides}{ephemeridesType},
and \configFile{satelliteModel}{satelliteModel} are only needed for some special parametrizations.
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

/** @brief Estimate accelerometer bias and scale.
* @ingroup programsGroup */
class InstrumentAccelerometerEstimateBiasScale
{
  InstrumentFile                 accFile;
  InstrumentFile                 accFileSim;
  InstrumentFile                 thrusterFile;
  InstrumentFile                 orbitFile;
  InstrumentFile                 starCameraFile;
  Double                         margin;
  EarthRotationPtr               earthRotation;
  EphemeridesPtr                 ephemerides;
  ParametrizationAccelerationPtr parameterAcceleration;
  SatelliteModelPtr              satellite;

  void observationEquation(UInt arcNo, AccelerometerArc &acc, Vector &l, Matrix &A, Matrix &B);

public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(InstrumentAccelerometerEstimateBiasScale, PARALLEL, "estimate accelerometer bias and scale.", Instrument)

/***********************************************/

void InstrumentAccelerometerEstimateBiasScale::run(Config &config)
{
  try
  {
    FileName accelerometerOutName, solutionOutName;
    FileName accelerometerInName,  accelerometerInNameSim;
    FileName thrusterInName;
    FileName orbitInName;
    FileName starCameraInName;
    FileName satelliteModelInName;

    renameDeprecatedConfig(config, "satelliteModel", "inputfileSatelliteModel",     date2time(2020, 8, 19));
    renameDeprecatedConfig(config, "parameter",      "parametrizationAcceleration", date2time(2020, 6,  3));

    readConfig(config, "outputfileAccelerometer",     accelerometerOutName,   Config::MUSTSET,  "",     "");
    readConfig(config, "outputfileSolution",          solutionOutName,        Config::OPTIONAL, "",     "");
    readConfig(config, "inputfileAccelerometer",      accelerometerInName,    Config::MUSTSET,  "",     "");
    readConfig(config, "inputfileAccelerometerSim",   accelerometerInNameSim, Config::MUSTSET,  "",     "");
    readConfig(config, "inputfileThruster",           thrusterInName,         Config::OPTIONAL, "",     "remove thruster events");
    readConfig(config, "marginThruster",              margin,                 Config::DEFAULT,  "1e-5", "margin size (on both sides) [seconds]");
    readConfig(config, "inputfileOrbit",              orbitInName,            Config::OPTIONAL, "",     "");
    readConfig(config, "inputfileStarCamera",         starCameraInName,       Config::OPTIONAL, "",     "");
    readConfig(config, "earthRotation",               earthRotation,          Config::OPTIONAL, "",     "");
    readConfig(config, "ephemerides",                 ephemerides,            Config::OPTIONAL, "jpl",  "may be needed by parametrizationAcceleration");
    readConfig(config, "inputfileSatelliteModel",     satelliteModelInName,   Config::OPTIONAL, "{groopsDataDir}/satelliteModel/", "satellite macro model");
    readConfig(config, "parametrizationAcceleration", parameterAcceleration,  Config::MUSTSET,  "",     "");
    if(isCreateSchema(config)) return;

     // ======================================================

    // init instrument files + satellite model
    // ---------------------------------------
    logStatus<<"read instrument data <"<<accelerometerInName<<">"<<Log::endl;
    logStatus<<"read simulated instrument data <"<<accelerometerInNameSim<<">"<<Log::endl;
    accFile.open(accelerometerInName);
    accFileSim.open(accelerometerInNameSim);
    orbitFile.open(orbitInName);
    starCameraFile.open(starCameraInName);
    thrusterFile.open(thrusterInName);
    InstrumentFile::checkArcCount({accFile, accFileSim, orbitFile, starCameraFile});

    if(!satelliteModelInName.empty())
      readFileSatelliteModel(satelliteModelInName, satellite);
    const UInt countParameter = parameterAcceleration->parameterCount();

     // ======================================================

    // estiamte accelerometer bias/scale
    // ---------------------------------
    Vector x;
    if(countParameter)
    {
      logStatus<<"estimate bias/scale"<<Log::endl;
      logInfo<<"  parameter count = "<<countParameter<<Log::endl;
      Matrix  N(countParameter, Matrix::SYMMETRIC);
      Vector  n(countParameter);
      Parallel::forEach(accFile.arcCount(), [&](UInt arcNo)
      {
        AccelerometerArc acc;
        Vector           l;
        Matrix           A, B;
        observationEquation(arcNo, acc, l, A, B);

        // remove thruster events (incl. margin)
        Arc thruster = thrusterFile.readArc(arcNo);
        if(thruster.size())
        {
          UInt idxThruster = 0;
          for(UInt i=0; i<acc.size(); i++)
          {
            while(idxThruster<thruster.size() && ((thruster.at(idxThruster).time-acc.at(i).time).seconds() < -margin))
              idxThruster++;
            if(idxThruster<thruster.size() && ((thruster.at(idxThruster).time-acc.at(i).time).seconds() > margin))
            {
              if(l.size()) l.row(3*i,3).setNull();
              if(A.size()) A.row(3*i,3).setNull();
              if(B.size()) B.row(3*i,3).setNull();
            }
          }
        }

        // accumulate normals
        // ------------------
        if(B.size())
          eliminationParameter(B, A, l);
        if(A.size())
        {
          rankKUpdate(1, A, N);
          matMult(1., A.trans(), l, n);
        }
      });
      Parallel::reduceSum(N);
      Parallel::reduceSum(n);

      if(Parallel::isMaster())
      {
        // regularize not used parameters
        // ------------------------------
        for(UInt i=0; i<N.rows(); i++)
          if(N(i,i)==0)
            N(i,i) = 1.0;

        x = solve(N, n);

        if(!solutionOutName.empty())
        {
          logStatus<<"write solution to <"<<solutionOutName<<">"<<Log::endl;
          writeFileMatrix(solutionOutName, x);
        }
      } // if(Parallel::isMaster())
      Parallel::broadCast(x);
    }

    // Apply estimated parameters
    // --------------------------
    std::vector<Arc> arcList(accFile.arcCount());
    Parallel::forEach(arcList, [&](UInt arcNo)
    {
      AccelerometerArc acc;
      Vector           l;
      Matrix           A, B;
      observationEquation(arcNo, acc, l, A, B);

      Vector Ax(l.rows());
      if(A.size())
        Ax = A*x;
      if(B.size())
        matMult(1, B, leastSquares(Matrix(B), l-Ax), Ax);

      for(UInt i=0; i<acc.size(); i++)
      {
        acc.at(i).acceleration.x() += Ax(3*i+0);
        acc.at(i).acceleration.y() += Ax(3*i+1);
        acc.at(i).acceleration.z() += Ax(3*i+2);
      }

      return acc;
    });

    if(Parallel::isMaster())
    {
      logStatus<<"write calibrated accelerometer file <"<<accelerometerOutName<<">"<<Log::endl;
      InstrumentFile::write(accelerometerOutName, arcList);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void InstrumentAccelerometerEstimateBiasScale::observationEquation(UInt arcNo, AccelerometerArc &acc, Vector &l, Matrix &A, Matrix &B)
{
  try
  {
    acc                         = accFile.readArc(arcNo);
    AccelerometerArc accSim     = accFileSim.readArc(arcNo);
    OrbitArc         orbit      = orbitFile.readArc(arcNo);
    StarCameraArc    starCamera = starCameraFile.readArc(arcNo);
    const UInt       epochCount = acc.size();
    Arc::checkSynchronized({acc, accSim, orbit, starCamera});

    parameterAcceleration->setIntervalArc(acc.at(0).time, acc.back().time+medianSampling(acc.times()));

    l = Vector(3*epochCount);
    A = Matrix(3*epochCount, parameterAcceleration->parameterCount());
    B = Matrix(3*epochCount, parameterAcceleration->parameterCountArc());

    for(UInt i=0; i<epochCount; i++)
    {
      l(3*i+0) = accSim.at(i).acceleration.x()-acc.at(i).acceleration.x();
      l(3*i+1) = accSim.at(i).acceleration.y()-acc.at(i).acceleration.y();
      l(3*i+2) = accSim.at(i).acceleration.z()-acc.at(i).acceleration.z();

      Rotary3d rotEarth, rotSat;
      Vector3d position, velocity;
      if(earthRotation)     rotEarth = earthRotation->rotaryMatrix(acc.at(i).time);
      if(starCamera.size()) rotSat   = starCamera.at(i).rotary;
      if(orbit.size())      position = orbit.at(i).position;
      if(orbit.size())      velocity = orbit.at(i).velocity;

      Matrix AEpoch(3, parameterAcceleration->parameterCount());
      Matrix BEpoch(3, parameterAcceleration->parameterCountArc());
      parameterAcceleration->compute(satellite, acc.at(i).time, position, velocity, rotSat, rotEarth, ephemerides, AEpoch, BEpoch);

      const Matrix R = (rotEarth*rotSat).matrix();  // rotate into accelerometer frame
      if(A.size()) matMult(1.0, R.trans(), AEpoch, A.row(3*i,3));
      if(B.size()) matMult(1.0, R.trans(), BEpoch, B.row(3*i,3));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
