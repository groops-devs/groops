/***********************************************/
/**
* @file graceAntennaCenterCorrectionArcCovariance.cpp
*
* @brief Compute antenna center correction from orbit configuration.
*
* @author Torsten Mayer-Guerr
* @date 2018-05-11
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program computes covariance information for the non­-stationary noise of the KBR antenna offset correction (AOC)
from the orientation covariance matrices provided in Level-1B products via variance propagation.
By using the output \configFile{outputfileSatelliteTrackingCovariance}{matrix} in \program{PreprocessingSst},
noise model distinguishes between the stationary noise of ranging observations and the non­stationary AOC noise.

The covariances are derived from the partial derivative of the AOC w.r.t. the roll/pitch/yaw rotations
and star camera covariances \configFile{inputfileScaCovariance1}{matrix} and \configFile{inputfileScaCovariance2}{matrix}.

The covariances for the range-rates and range-acceleration are computed by differentiating
an interpolation polynomial of degree \config{interpolationDegree}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/polynomial.h"
#include "parallel/matrixDistributed.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"
#include "misc/grace/graceKBandGeometry.h"

/***** CLASS ***********************************/

/** @brief Compute antenna center correction from orbit configuration.
* @ingroup programsGroup */
class GraceAntennaCenterCorrectionArcCovariance
{
  void rotaryCholesky(const std::vector<Time> &times, const Covariance3dArc &starCameraCovariance, const Vector &sigmaAxisAcc, UInt degree, MatrixDistributed &normals) const;

public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(GraceAntennaCenterCorrectionArcCovariance, PARALLEL, "compute antenna center correction from orbit configuration", Grace, Covariance)

/***********************************************/

void GraceAntennaCenterCorrectionArcCovariance::run(Config &config)
{
  try
  {
    FileName    fileNameCovariance;
    FileName    fileNameOrbit1, fileNameOrbit2;
    FileName    fileNameStarCamera1, fileNameStarCamera2;
    FileName    fileNameScaCov1, fileNameScaCov2;
    UInt        sstType;
    Vector      sigmaAxisAcc(3);
    Vector3d    center1, center2;
    UInt        degree;
    std::string choice;

    readConfig(config, "outputfileSatelliteTrackingCovariance", fileNameCovariance, Config::OPTIONAL, "", "corrections for range, range-rate, and range-accelerations");
    if(readConfigChoice(config, "sstType", choice, Config::MUSTSET, "", ""))
    {
      if(readConfigChoiceElement(config, "range",             choice, "")) sstType = 0;
      if(readConfigChoiceElement(config, "rangeRate",         choice, "")) sstType = 1;
      if(readConfigChoiceElement(config, "rangeAcceleration", choice, "")) sstType = 2;
      endChoice(config);
    }
    readConfig(config, "inputfileOrbit1",         fileNameOrbit1,      Config::MUSTSET,  "", "");
    readConfig(config, "inputfileOrbit2",         fileNameOrbit2,      Config::MUSTSET,  "", "");
    readConfig(config, "inputfileStarCamera1",    fileNameStarCamera1, Config::MUSTSET,  "", "");
    readConfig(config, "inputfileStarCamera2",    fileNameStarCamera2, Config::MUSTSET,  "", "");
    readConfig(config, "inputfileScaCovariance1", fileNameScaCov1,     Config::MUSTSET,  "", "");
    readConfig(config, "inputfileScaCovariance2", fileNameScaCov2,     Config::MUSTSET,  "", "");
    readConfig(config, "sigmaAccelerometerX",     sigmaAxisAcc(0),     Config::DEFAULT,  "10e-7", "[rad/s^2]");
    readConfig(config, "sigmaAccelerometerY",     sigmaAxisAcc(1),     Config::DEFAULT,  "2e-7",  "[rad/s^2]");
    readConfig(config, "sigmaAccelerometerZ",     sigmaAxisAcc(2),     Config::DEFAULT,  "2e-7",  "[rad/s^2]");
    if(readConfigChoice(config, "antennaCenters",     choice,          Config::MUSTSET,  "",      "KBR antenna phase center"))
    {
      if(readConfigChoiceElement(config, "value", choice, ""))
      {
        readConfig(config, "center1X", center1.x(), Config::DEFAULT,   "1.4451172588", "x-coordinate of antenna position in SRF [m] for GRACEA");
        readConfig(config, "center1Y", center1.y(), Config::DEFAULT,  "-0.0004233040", "y-coordinate of antenna position in SRF [m] for GRACEA");
        readConfig(config, "center1Z", center1.z(), Config::DEFAULT,   "0.0022786600", "z-coordinate of antenna position in SRF [m] for GRACEA");
        readConfig(config, "center2X", center2.x(), Config::DEFAULT,   "1.4443870350", "x-coordinate of antenna position in SRF [m] for GRACEB");
        readConfig(config, "center2Y", center2.y(), Config::DEFAULT,   "0.0005761203", "y-coordinate of antenna position in SRF [m] for GRACEB");
        readConfig(config, "center2Z", center2.z(), Config::DEFAULT,   "0.0033040887", "z-coordinate of antenna position in SRF [m] for GRACEB");
      }
      if(readConfigChoiceElement(config, "file",  choice, ""))
      {
        FileName fileName;
        readConfig(config, "inputAntennaCenters", fileName, Config::MUSTSET, "", "");
        if(!isCreateSchema(config))
        {
          Matrix x;
          readFileMatrix(fileName, x);
          center1 = Vector3d(x(0,0), x(1,0), x(2,0));
          center2 = Vector3d(x(3,0), x(4,0), x(5,0));
        }
      }
      endChoice(config);
    }
    readConfig(config, "interpolationDegree", degree, Config::DEFAULT,  "2", "differentiation by polynomial approximation of degree n");
    if(isCreateSchema(config)) return;

    // =======================================================

    logStatus<<"read orbit and star camera data and generate antenna offset corrections"<<Log::endl;
    InstrumentFile  orbit1File(fileNameOrbit1);
    InstrumentFile  orbit2File(fileNameOrbit2);
    InstrumentFile  starCamera1File(fileNameStarCamera1);
    InstrumentFile  starCamera2File(fileNameStarCamera2);
    InstrumentFile  starCameraCovariance1File(fileNameScaCov1);
    InstrumentFile  starCameraCovariance2File(fileNameScaCov2);
    InstrumentFile::checkArcCount({orbit1File, orbit2File, starCamera1File, starCamera2File});

    Parallel::forEach(orbit1File.arcCount(), [&](UInt arcNo)
    {
      OrbitArc        orbit1                = orbit1File.readArc(arcNo);
      OrbitArc        orbit2                = orbit2File.readArc(arcNo);
      StarCameraArc   starCamera1           = starCamera1File.readArc(arcNo);
      StarCameraArc   starCamera2           = starCamera2File.readArc(arcNo);
      Covariance3dArc starCameraCovariance1 = starCameraCovariance1File.readArc(arcNo);
      Covariance3dArc starCameraCovariance2 = starCameraCovariance2File.readArc(arcNo);
      Arc::checkSynchronized({orbit1, orbit2, starCamera1, starCamera2});

      const std::vector<Time> times = orbit1.times();
      const UInt   epochCount = times.size();
      const Double dt = medianSampling(times).seconds();

      MatrixDistributed normals[2];
      rotaryCholesky(times, starCameraCovariance1, sigmaAxisAcc, degree, normals[0]);
      rotaryCholesky(times, starCameraCovariance2, sigmaAxisAcc, degree, normals[1]);

      Matrix SparseJacobian[2];
      GraceKBandGeometry::partialOfAntennaCenterCorrectionWrtRollPitchYaw(orbit1, orbit2, starCamera1, starCamera2, center1, center2, SparseJacobian[0], SparseJacobian[1]);

      Matrix Cov(epochCount, Matrix::SYMMETRIC);
      Polynomial p(degree);
      for(UInt id=0; id<2; id++) // GRACE A/B
      {
        Matrix A(epochCount, 3*epochCount);
        for(UInt i=0; i<epochCount; i++)
          copy(SparseJacobian[id].row(i), A.slice(i, 3*i, 1, 3));

        if(sstType == 0)      // range
          A = A.trans();
        else if(sstType == 1) // range rate
          A = p.derivative(dt, A).trans();
        else if(sstType == 2) // range acceleration
          A = p.derivative2nd(dt, A).trans();

        // apply covariance
        normals[id].triangularTransSolve(A);
        rankKUpdate(1., A, Cov);
      }

      writeFileMatrix(fileNameCovariance.appendBaseName(".arc"+arcNo%"%03i"s),  Cov);
    }); // forEach
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GraceAntennaCenterCorrectionArcCovariance::rotaryCholesky(const std::vector<Time> &times, const Covariance3dArc &starCameraCovariance, const Vector &sigmaAxisAcc, UInt degree, MatrixDistributed &normals) const
{
  try
  {
    const UInt   epochCount = times.size();
    const Double dt         = medianSampling(times).seconds();

    std::vector<UInt> blockIndex(1, 0);
    for(UInt i=0; i<epochCount; i++)
      blockIndex.push_back( blockIndex.back()+3 );
    normals.initEmpty(blockIndex, Parallel::selfCommunicator());
    for(UInt i=0; i<epochCount; i++)
      normals.setBlock(i, i);

    // polynomial interpolation matrix
    // -------------------------------
    Matrix P(degree+1, degree+1);
    for(UInt i=0; i<P.rows(); i++)
    {
      P(0,i) = 1.0;
      for(UInt n=1; n<P.columns(); n++)
        P(n,i) = (i-degree/2.) * P(n-1,i);
    }
    inverse(P);

    // star camera observations
    // ------------------------
    UInt idxSca = 0;
    const std::vector<Time> timesCov = starCameraCovariance.times();
    for(UInt i=0; i<epochCount; i++)
    {
      // search time slot
      while((idxSca < timesCov.size()) && (timesCov.at(idxSca) < times.at(i)))
        idxSca++;
      if(idxSca >= timesCov.size())
        break;
      if(timesCov.at(idxSca) > times.at(i))
        continue;

      normals.N(i,i) = starCameraCovariance.at(idxSca).covariance.matrix();
      inverse(normals.N(i,i));
    } // for(i=epochCount)

    // angular accelerations observations
    // ----------------------------------
    for(UInt i=0; i<epochCount; i++)
    {
      // polynomial 2nd differentiation
      const UInt idx = std::min(std::max(i, degree/2)-degree/2, epochCount-degree-1); // interval
      Vector factors = Vector(degree+1);  // derivation coefficients
      const Double tau = i-idx-degree/2.;
      Double f = 1./(dt*dt);
      for(UInt n=2; n<factors.rows(); n++)
      {
        axpy(n*(n-1)*f, P.column(n), factors);
        f *= tau;
      }

      // dacc(i)/dq(k)
      std::vector<Matrix> accA(factors.size());
      for(UInt k=0; k<factors.size(); k++) // w_dot(i) = sum_k factor(k) alpha(k)
        accA.at(k) = factors(k) * identityMatrix(3);

      // decorrelation
      for(UInt idAxis=0; idAxis<3; idAxis++)
        for(UInt k=0; k<accA.size(); k++)
          accA.at(k).row(idAxis) *= 1./sigmaAxisAcc(idAxis);

      // accumulate normals
      for(UInt k=0; k<factors.size(); k++)
      {
        rankKUpdate(1., accA.at(k), normals.N(idx+k, idx+k)); // N11(i)
        for(UInt j=k+1; j<factors.size(); j++)                // N11(i), secondary diagonal
        {
          normals.setBlock(idx+k, idx+j);
          matMult(1., accA.at(k).trans(), accA.at(j), normals.N(idx+k, idx+j));
        }
      }
    } // for(i=epochCount)

    // cholesky
    // --------
    normals.cholesky(FALSE/*timing*/);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
