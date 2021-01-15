/***********************************************/
/**
* @file instrumentStarCameraAngularAccelerometerFusion.cpp
*
* @brief Estimation of the satellites orientation from star camera and angular accelerometer data.
*
* @author Beate Klinger
* @author Torsten Mayer-Guerr
* @date 2017-06-22
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program estimates the satellites orientation from star camera data
\configFile{inputfileStarCamera}{instrument} and angular accelerometer data
\configFile{inputfileAngularAcc}{instrument}. The combination of both observation types
is achieved in a least square adjustment. The optimal weighting between the two different
observation groups is achieved by means of VCE in combination with a robust estimator.
The system of linearized observation equations within the sensor fusion approach can be formulated as:
\begin{equation}
  \begin{bmatrix}
  \M l_{ACC1B}\\
  \M l_{SCA1B}
  \end{bmatrix}
  =
  \begin{bmatrix}
  \M A_{ACC1B} & \M B_{ACC1B}\\
  \M A_{SCA1B} & \M 0
  \end{bmatrix}
  \begin{bmatrix}
  \M q\\
  \M b
  \end{bmatrix}
  =
  \begin{bmatrix}
  \frac{\partial \dot{\boldsymbol{\omega}}}{\partial \M q} & \frac{\partial \dot{\boldsymbol{\omega}}}{\partial \M b}\\
  \M I & \M 0
  \end{bmatrix}
  \begin{bmatrix}
  \M q\\
  \M b
  \end{bmatrix}
\end{equation}
with
\begin{equation}\begin{split}
  \M l_{ACC1B}  &= \dot{\boldsymbol{\omega}}_{ACC1B} - \dot{\boldsymbol{\omega}}_{0}, \\
  \M l_{SCA1B}  &= \M q_{SCA1B} - \M q_{0}, \\
  \M q_{Fusion} &= \M q + \M q_{0}.
\end{split}\end{equation}
The reference values $\M q_{0}$ and $\dot{\boldsymbol{\omega}}_{0}$ are derived
from \configFile{inputfileStarCameraReference}{instrument}. In the course of the estimation,
the accelerometer data is calibrated, by setting a bias factor $\M b$ with \config{accBias}.
)";

/***********************************************/

#include "programs/program.h"
#include "parallel/matrixDistributed.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"
#include "misc/varianceComponentEstimation.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"

/***** CLASS ***********************************/

/** @brief Estimation of the satellites orientation from star camera and angular accelerometer data.
* @ingroup programsGroup */
class InstrumentStarCameraAngularAccelerometerFusion
{
public:

  class ArcObservation
  {
  public:
    StarCameraArc    starCamera;
    Covariance3dArc  starCameraCovariance;
    AccelerometerArc angularAcc;
    MiscValueArc     epochSigmaSca, epochSigmaAcc;
    Vector           sigmaSca, sigmaAcc;
    Vector           x;
  };

  FileName                   fileNameOutStarCamera;
  FileName                   fileNameOutCovariance;
  FileName                   fileNameOutCovarianceMatrix;
  FileName                   fileNameOutAngularAcc;
  FileName                   fileNameOutEpochSigmasStarCamera;
  FileName                   fileNameOutEpochSigmasAccelerometer;
  InstrumentFile             starCamera0File;
  std::vector<Time>          timesSca;
  Matrix                     quaternionSca;
  Covariance3dArc            starCameraCovariance;
  std::vector<Time>          timesAcc;
  AccelerometerArc           angularAcc;
  Bool                       correctAccNonSquare;
  ParametrizationTemporalPtr accBias, accScale;
  Bool                       perAxisSca, perAxisAcc;
  Double                     sigmaSca0;
  Double                     sigmaAccX, sigmaAccY, sigmaAccZ;
  Double                     huber, huberPower;
  UInt                       degree;
  Matrix                     P; // polynomial interpolation
  UInt                       iterCount;


  Bool           index(const Time &time, std::vector<Time> &times, UInt &idx) const;
  ArcObservation computeArc(UInt arc);
  static Vector  computeRedundancy(std::vector<Matrix> &A, MatrixDistributed &normals, Double threshold=1e-4);
  static void    estimateAccuracy(std::vector<Vector> &e, std::vector<Vector> &redundancy, Bool sigmaPerAxis, Double huber, Double huberPower, Vector &sigmaAxis, Vector &sigmaEpoch);

  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(InstrumentStarCameraAngularAccelerometerFusion, PARALLEL, "estimation of the satellites orientation from star camera and angular accelerometer data.", Grace, Instrument)

/***********************************************/

template<> void save(OutArchive &oa, const InstrumentStarCameraAngularAccelerometerFusion::ArcObservation &arc)
{
  oa<<nameValue("starCamera",           arc.starCamera);
  oa<<nameValue("starCameraCovariance", arc.starCameraCovariance);
  oa<<nameValue("angularAcc",           arc.angularAcc);
  oa<<nameValue("sigmaSca",             arc.sigmaSca);
  oa<<nameValue("sigmaAcc",             arc.sigmaAcc);
  oa<<nameValue("epochSigmaSca",        arc.epochSigmaSca);
  oa<<nameValue("epochSigmaAcc",        arc.epochSigmaAcc);
  oa<<nameValue("x",                    arc.x);
}

template<> void load(InArchive  &ia, InstrumentStarCameraAngularAccelerometerFusion::ArcObservation &arc)
{
  ia>>nameValue("starCamera",           arc.starCamera);
  ia>>nameValue("starCameraCovariance", arc.starCameraCovariance);
  ia>>nameValue("angularAcc",           arc.angularAcc);
  ia>>nameValue("sigmaSca",             arc.sigmaSca);
  ia>>nameValue("sigmaAcc",             arc.sigmaAcc);
  ia>>nameValue("epochSigmaSca",        arc.epochSigmaSca);
  ia>>nameValue("epochSigmaAcc",        arc.epochSigmaAcc);
  ia>>nameValue("x",                    arc.x);
}

/***********************************************/

void InstrumentStarCameraAngularAccelerometerFusion::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName fileNameOutSolution;
    FileName fileNameStarCamera0;
    FileName fileNameStarCamera, fileNameStarCameraCov, fileNameAngularAcc;

    readConfig(config, "outputfileStarCamera",              fileNameOutStarCamera,               Config::MUSTSET,  "",      "combined quaternions");
    readConfig(config, "outputfileCovariance",              fileNameOutCovariance,               Config::OPTIONAL, "",      "epoch-wise covariance matrix");
    readConfig(config, "outputfileCovarianceMatrix",        fileNameOutCovarianceMatrix,         Config::OPTIONAL, "",      "full arc-wise covariance matrix per arc. arc number is appended to filename");
    readConfig(config, "outputfileEpochSigmaStarCamera",    fileNameOutEpochSigmasStarCamera,    Config::OPTIONAL, "",      "from vce and outlier detection");
    readConfig(config, "outputfileEpochSigmaAccelerometer", fileNameOutEpochSigmasAccelerometer, Config::OPTIONAL, "",      "from vce and outlier detection");
    readConfig(config, "outputfileAngularAcc",              fileNameOutAngularAcc,               Config::OPTIONAL, "",      "angular acceleration observations (bias removed)");
    readConfig(config, "outputfileSolution",                fileNameOutSolution,                 Config::OPTIONAL, "",      "estimated parameter (one column for each arc)");
    readConfig(config, "inputfileStarCameraReference",      fileNameStarCamera0,                 Config::MUSTSET,  "",      "quaternions as taylor point");
    readConfig(config, "inputfileStarCamera",               fileNameStarCamera,                  Config::MUSTSET,  "",      "star camera observations");
    readConfig(config, "inputfileStarCameraCovariance",     fileNameStarCameraCov,               Config::OPTIONAL, "",      "star camera observations");
    readConfig(config, "inputfileAngularAcc",               fileNameAngularAcc,                  Config::MUSTSET,  "",      "angular acceleration observations");
    readConfig(config, "correctAccNonQuadratic",            correctAccNonSquare,                 Config::DEFAULT,  "1",     "apply correction (non-square proof mass)");
    readConfig(config, "accBias",                           accBias,                             Config::DEFAULT,  "",      "accelerometer bias per interval and axis");
    readConfig(config, "accScale",                          accScale,                            Config::DEFAULT,  "",      "accelerometer scale per interval and axis");
    readConfig(config, "sigmaStarcamera",                   sigmaSca0,                           Config::DEFAULT,  "1.",    "[rad]");
    readConfig(config, "sigmaAccelerometerX",               sigmaAccX,                           Config::DEFAULT,  "10e-7", "[rad/s^2]");
    readConfig(config, "sigmaAccelerometerY",               sigmaAccY,                           Config::DEFAULT,  "2e-7",  "[rad/s^2]");
    readConfig(config, "sigmaAccelerometerZ",               sigmaAccZ,                           Config::DEFAULT,  "2e-7",  "[rad/s^2]");
    readConfig(config, "estimateSigmaScaPerAxis",           perAxisSca,                          Config::DEFAULT,  "0",     "separate variance factor for roll, pitch, yaw, instead of one common factor.""");
    readConfig(config, "estimateSigmaAccPerAxis",           perAxisAcc,                          Config::DEFAULT,  "1",     "separate variance factor for each accelerometer axis, instead of one common factor.");
    readConfig(config, "huber",                             huber,                               Config::DEFAULT,  "2.5",   "residuals > huber*sigma0 are downweighted");
    readConfig(config, "huberPower",                        huberPower,                          Config::DEFAULT,  "1.5",   "residuals > huber: sigma=(e/huber)^power*sigma0");
    readConfig(config, "interpolationDegree",               degree,                              Config::DEFAULT,  "2",     "");
    readConfig(config, "iterationCount",                    iterCount,                           Config::DEFAULT,  "15",    "non linear equation solved iteratively");
    if(isCreateSchema(config)) return;

    // Read instruments
    // ----------------
    StarCameraArc starCamera = InstrumentFile::read(fileNameStarCamera);
    starCameraCovariance     = InstrumentFile::read(fileNameStarCameraCov);
    Arc::checkSynchronized({starCamera, starCameraCovariance});
    angularAcc               = InstrumentFile::read(fileNameAngularAcc);
    timesSca                 = starCamera.times();
    quaternionSca            = starCamera.matrix().column(1,4); // without time
    timesAcc                 = angularAcc.times();

    // polynomial interpolation matrix
    // -------------------------------
    P = Matrix(degree+1, degree+1);
    for(UInt i=0; i<P.rows(); i++)
    {
      P(0,i) = 1.0;
      for(UInt n=1; n<P.columns(); n++)
        P(n,i) = (i-degree/2.) * P(n-1,i);
    }
    inverse(P);

    // computing quaternions
    // ---------------------
    logStatus<<"computing quaternions"<<Log::endl;
    starCamera0File.open(fileNameStarCamera0);
    std::vector<ArcObservation> arcs(starCamera0File.arcCount());
    Parallel::forEach(arcs, [this] (UInt arcNo) {return computeArc(arcNo);}, comm);

    if(Parallel::isMaster(comm))
    {
      for(UInt arcNo=0; arcNo<arcs.size(); arcNo++)
        logInfo<<"  "<<arcNo%"%2i. arc "s
               <<"sigmaScaX = "<<arcs.at(arcNo).sigmaSca(0)%"%8.2e, "s
               <<"sigmaScaY = "<<arcs.at(arcNo).sigmaSca(1)%"%8.2e, "s
               <<"sigmaScaZ = "<<arcs.at(arcNo).sigmaSca(2)%"%8.2e, "s
               <<"sigmaAccX = "<<arcs.at(arcNo).sigmaAcc(0)%"%8.2e, "s
               <<"sigmaAccY = "<<arcs.at(arcNo).sigmaAcc(1)%"%8.2e, "s
               <<"sigmaAccZ = "<<arcs.at(arcNo).sigmaAcc(2)%"%8.2e"s<<Log::endl;
    }

    if(Parallel::isMaster(comm) && !fileNameOutStarCamera.empty())
    {
      logStatus<<"write starCamera file <"<<fileNameOutStarCamera<<">"<<Log::endl;
      std::vector<Arc> arcList(arcs.size());
      for(UInt arcNo=0; arcNo<arcs.size(); arcNo++)
        arcList.at(arcNo) = arcs.at(arcNo).starCamera;
      InstrumentFile::write(fileNameOutStarCamera, arcList);
    }

    if(Parallel::isMaster(comm) && !fileNameOutCovariance.empty())
    {
      logStatus<<"write covariance file <"<<fileNameOutCovariance<<">"<<Log::endl;
      std::vector<Arc> arcList(arcs.size());
      for(UInt arcNo=0; arcNo<arcs.size(); arcNo++)
        arcList.at(arcNo) = arcs.at(arcNo).starCameraCovariance;
      InstrumentFile::write(fileNameOutCovariance, arcList);
    }

    if(Parallel::isMaster(comm) && !fileNameOutAngularAcc.empty())
    {
      logStatus<<"write angular acceleration file <"<<fileNameOutAngularAcc<<">"<<Log::endl;
      std::vector<Arc> arcList(arcs.size());
      for(UInt arcNo=0; arcNo<arcs.size(); arcNo++)
        arcList.at(arcNo) = arcs.at(arcNo).angularAcc;
      InstrumentFile::write(fileNameOutAngularAcc, arcList);
    }

    if(Parallel::isMaster(comm) && !fileNameOutSolution.empty())
    {
      logStatus<<"write solution file <"<<fileNameOutSolution<<">"<<Log::endl;
      Matrix x(arcs.at(0).x.rows(), arcs.size());
      for(UInt arcNo=0; arcNo<arcs.size(); arcNo++)
        copy(arcs.at(arcNo).x, x.column(arcNo));
      writeFileMatrix(fileNameOutSolution, x);
    }

    if(Parallel::isMaster(comm) && !fileNameOutEpochSigmasStarCamera.empty())
    {
      logStatus<<"write star camera epoch sigmas to file <"<<fileNameOutEpochSigmasStarCamera<<">"<<Log::endl;
      std::vector<Arc> arcList(arcs.size());
      for(UInt arcNo=0; arcNo<arcs.size(); arcNo++)
        arcList.at(arcNo) = arcs.at(arcNo).epochSigmaSca;
      InstrumentFile::write(fileNameOutEpochSigmasStarCamera, arcList);
    }

    if(Parallel::isMaster(comm) && !fileNameOutEpochSigmasAccelerometer.empty())
    {
      logStatus<<"write accelerometer epoch sigmas to file <"<<fileNameOutEpochSigmasAccelerometer<<">"<<Log::endl;
      std::vector<Arc> arcList(arcs.size());
      for(UInt arcNo=0; arcNo<arcs.size(); arcNo++)
        arcList.at(arcNo) = arcs.at(arcNo).epochSigmaAcc;
      InstrumentFile::write(fileNameOutEpochSigmasAccelerometer, arcList);
    }

  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool InstrumentStarCameraAngularAccelerometerFusion::index(const Time &time, std::vector<Time> &times, UInt &idx) const
{
  try
  {
    const Time timeMargin = seconds2time(1e-5);
    // search time slot
    while((idx < times.size()) && (times.at(idx)+timeMargin < time))
      idx++;
    if(idx >= times.size())
      return FALSE;
    if(times.at(idx)-timeMargin > time)
      return FALSE;
    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

InstrumentStarCameraAngularAccelerometerFusion::ArcObservation InstrumentStarCameraAngularAccelerometerFusion::computeArc(UInt arcNo)
{
  try
  {
    const StarCameraArc     starCamera0 = starCamera0File.readArc(arcNo);
    const std::vector<Time> times       = starCamera0.times();
    const UInt              epochCount  = starCamera0.size();
    const Double            dt          = medianSampling(times).seconds();

    // approximate values
    // ------------------
    Matrix q = starCamera0.matrix().column(1,4);
    // check sign
    UInt idxSca = 0;
    for(UInt i=0; i<epochCount; i++)
      if(index(times.at(i), timesSca, idxSca))
        if(inner(quaternionSca.row(idxSca), q.row(i)) < 0)
        {
          q *= -1;
          break;
        }

    // accelerometer calibration
    // -------------------------
    accBias->setInterval (times.at(0), times.back()+seconds2time(dt), TRUE/*estimatePerArc*/);
    accScale->setInterval(times.at(0), times.back()+seconds2time(dt), TRUE/*estimatePerArc*/);

    UInt countAccCal = 0;
    const UInt idxAccBias  = countAccCal; countAccCal += 3 * accBias->parameterCount();
    const UInt idxAccScale = countAccCal; countAccCal += 3 * accScale->parameterCount();

    // epoch wise sigmas
    // -----------------
    Vector sigmaAxisSca = {sigmaSca0, sigmaSca0, sigmaSca0};
    Vector sigmaEpochSca(epochCount);
    for(UInt i=0; i<epochCount; i++)
      sigmaEpochSca(i) = 1.;

    Vector sigmaAxisAcc = {sigmaAccX, sigmaAccY, sigmaAccZ};
    Vector sigmaEpochAcc(epochCount);
    for(UInt i=0; i<epochCount; i++)
      sigmaEpochAcc(i) = 1.;

    // normal equations (kite structure)
    // ---------------------------------
    Double lPl      = 0;
    UInt   obsCount = 0;
    Vector n(3*epochCount+countAccCal);
    Vector xAccCal(countAccCal);

    std::vector<UInt> blockIndex(1, 0);
    for(UInt i=0; i<epochCount; i++)
      blockIndex.push_back( blockIndex.back()+3 );

    const UInt idxBlockCal = epochCount;
    if(countAccCal)
      blockIndex.push_back( blockIndex.back()+countAccCal );

    MatrixDistributed normals;
    normals.initEmpty(blockIndex, Parallel::selfCommunicator());
    const UInt blockCount = normals.blockCount();
    for(UInt i=0; i<epochCount; i++)
      normals.setBlock(i, i);

    // =============================================

    Double sigma0 = 1.;
    std::vector<Bool> hasSca(epochCount, FALSE);
    std::vector<Bool> hasAcc(epochCount, FALSE);
    for(UInt iter=0; iter<iterCount; iter++)
    {
      // init normal equations
      // ---------------------
      normals.setNull();
      n.setNull();
      lPl      = 0;
      obsCount = 0;

      std::vector<Matrix> Q(epochCount);
      for(UInt i=0; i<epochCount; i++)
        Q.at(i) = {{-q(i,1),  q(i,0),  q(i,3), -q(i,2)},
                   {-q(i,2), -q(i,3),  q(i,0),  q(i,1)},
                   {-q(i,3),  q(i,2), -q(i,1),  q(i,0)}};

      // star camera observations
      // ------------------------
      std::vector<Matrix> scaA(epochCount);
      std::vector<Vector> scaL(epochCount);
      UInt                idxSca = 0;
      for(UInt i=0; i<epochCount; i++)
        if(index(times.at(i), timesSca, idxSca))
        {
          hasSca.at(i) = TRUE;

          scaA.at(i) = identityMatrix(3);
          scaL.at(i) = 2*Q.at(i)*quaternionSca.row(idxSca).trans();

          // decorrelation
          Matrix W = (starCameraCovariance.size()) ? starCameraCovariance.at(idxSca).covariance.matrix() : identityMatrix(3, Matrix::SYMMETRIC);
          cholesky(W);
          for(UInt idAxis=0; idAxis<3; idAxis++)
            W.column(idAxis) *= sigmaEpochSca(i)*sigmaAxisSca(idAxis);
          triangularSolve(1., W.trans(), scaA.at(i));
          triangularSolve(1., W.trans(), scaL.at(i));

          matMult(1., scaA.at(i).trans(), scaL.at(i), n.row(normals.blockIndex(i), normals.blockSize(i)));
          rankKUpdate(1., scaA.at(i), normals.N(i,i));
          lPl      += quadsum(scaL.at(i));
          obsCount += scaL.at(i).rows();
        } // for(i=epochCount)

      // angular accelerations observations
      // ----------------------------------
      std::vector<UInt>                accIdxA(epochCount); // start of A matrix
      std::vector<std::vector<Matrix>> accA(epochCount);    // quaternions
      std::vector<Matrix>              accB(epochCount);    // bias & scale
      std::vector<Vector>              accL(epochCount);
      UInt                             idxAcc = 0;
      for(UInt i=0; i<epochCount; i++)
        if(index(times.at(i), timesAcc, idxAcc))
        {
          hasAcc.at(i) = TRUE;

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

          // reduced observations
          accL.at(i) = angularAcc.at(idxAcc).acceleration.vector();
          for(UInt k=0; k<factors.size(); k++)
            matMult(-2*factors(k), Q.at(i), q.row(idx+k).trans(), accL.at(i));

          if(correctAccNonSquare)
          {
            // polynomial differentiation
            UInt idxRate = std::min(std::max(i, UInt(1))-1, epochCount-2);
            std::vector<Double> factorsRate = {-1/dt, 0, 1/dt};
            if((i == 0) || (i == epochCount-1))
              factorsRate = {-1/dt, 1/dt};

            Vector angularRate(3);
            for(UInt k=0; k<factorsRate.size(); k++)
              matMult(2*factorsRate.at(k), Q.at(i), q.row(idxRate+k).trans(), angularRate);

            const Double Ix = 1.037e-5;
            const Double Iy = 1.959e-5;
            const Double Iz = Ix;
            accL.at(i)(0) -= (Iz-Iy)/Ix * angularRate(2) * angularRate(1);
            accL.at(i)(1) -= (Iz-Ix)/Iy * angularRate(2) * angularRate(0);
            accL.at(i)(2) -= (Iy-Ix)/Iz * angularRate(1) * angularRate(0);
          }

          // dacc(i)/dq(k)
          accIdxA.at(i) = idx;
          accA.at(i).resize(factors.size());
          for(UInt k=0; k<factors.size(); k++) // w_dot(i) = sum_k factor(k) alpha(k)
            accA.at(i).at(k) = factors(k) * identityMatrix(3);

          // calibration parameter
          accB.at(i) = Matrix(3, countAccCal);
          if(accBias->parameterCount())
            accBias->designMatrix(times.at(i), identityMatrix(3), accB.at(i).column(idxAccBias, 3*accBias->parameterCount()));
          if(accScale->parameterCount())
          {
            Matrix S(3, 3);
            S(0, 0) = angularAcc.at(idxAcc).acceleration.x();
            S(1, 1) = angularAcc.at(idxAcc).acceleration.y();
            S(2, 2) = angularAcc.at(idxAcc).acceleration.z();
            accScale->designMatrix(times.at(i), S, accB.at(i).column(idxAccScale, S.columns()*accScale->parameterCount()));
          }

          // decorrelation
          for(UInt idAxis=0; idAxis<3; idAxis++)
          {
            accL.at(i)(idAxis) *= 1./(sigmaEpochAcc(i)*sigmaAxisAcc(idAxis));
            for(UInt k=0; k<accA.at(i).size(); k++)
              accA.at(i).at(k).row(idAxis) *= 1./(sigmaEpochAcc(i)*sigmaAxisAcc(idAxis));
            if(accB.at(i).size())
              accB.at(i).row(idAxis) *= 1./(sigmaEpochAcc(i)*sigmaAxisAcc(idAxis));
          }

          // accumulate normals
          for(UInt k=0; k<factors.size(); k++)
          {
            matMult(1., accA.at(i).at(k).trans(), accL.at(i), n.row(normals.blockIndex(idx+k), normals.blockSize(idx+k)));
            rankKUpdate(1., accA.at(i).at(k), normals.N(idx+k, idx+k));     // N11(i)
            for(UInt j=k+1; j<factors.size(); j++)                          // N11(i), secondary diagonal
            {
              normals.setBlock(idx+k, idx+j);
              matMult(1., accA.at(i).at(k).trans(), accA.at(i).at(j), normals.N(idx+k, idx+j));
            }
            if(accB.at(i).size())
            {
              normals.setBlock(idx+k, idxBlockCal);
              matMult(1., accA.at(i).at(k).trans(), accB.at(i), normals.N(idx+k, idxBlockCal)); // N12(i)
            }
          }
          if(accB.at(i).size())
          {
            normals.setBlock(idxBlockCal, idxBlockCal);
            matMult(1., accB.at(i).trans(), accL.at(i), n.row(normals.blockIndex(idxBlockCal), normals.blockSize(idxBlockCal))); // n2
            rankKUpdate(1., accB.at(i), normals.N(idxBlockCal, idxBlockCal)); // N22
          }
          lPl += quadsum(accL.at(i));
          obsCount += accL.at(i).rows();
        } // for(i=epochCount)

      // regularize epochs without SCA
      // -----------------------------
      for(UInt i=0; i<epochCount; i++)
        if(!hasSca.at(i))
        {
          const Double sigma = 0.25*DEG2RAD;
          axpy(1./(sigma*sigma), identityMatrix(3), normals.N(i,i));
        }

      // Solve
      // -----
      const Vector dx = normals.solve(n, FALSE/*timing*/);
      const Double sigma0Old = sigma0;
      sigma0 = sqrt((lPl-inner(dx,n))/(obsCount-dx.rows()));
      if(countAccCal)
        xAccCal += dx.row(normals.blockIndex(idxBlockCal), normals.blockSize(idxBlockCal));

      // update approximate values
      // -------------------------
      for(UInt i=0; i<epochCount; i++)
      {
        matMult(0.5, Q.at(i).trans(), dx.row(normals.blockIndex(i), normals.blockSize(i)), q.row(i).trans());
        q.row(i) *= 1./norm(q.row(i));
      }

      // remove calibration parameter
      // ----------------------------
      if(countAccCal)
      {
        UInt idxAcc = 0;
        for(UInt i=0; i<epochCount; i++)
          if(index(times.at(i), timesAcc, idxAcc))
          {
            Vector accCal = accB.at(i) * dx.row(normals.blockIndex(idxBlockCal), normals.blockSize(idxBlockCal));
            angularAcc.at(idxAcc).acceleration.x() -= (sigmaEpochAcc(i)*sigmaAxisAcc(0)) * accCal(0);
            angularAcc.at(idxAcc).acceleration.y() -= (sigmaEpochAcc(i)*sigmaAxisAcc(1)) * accCal(1);
            angularAcc.at(idxAcc).acceleration.z() -= (sigmaEpochAcc(i)*sigmaAxisAcc(2)) * accCal(2);
          }
      }

      // Accuracy & outlier detection (star camera)
      // ------------------------------------------
      {
        std::vector<Vector> e(epochCount);
        std::vector<Vector> redundancy(epochCount);
        std::vector<Matrix> A(epochCount+1);
        for(UInt i=0; i<epochCount; i++)
          if(scaL.at(i).size())
          {
            // decorrelated residuals
            e.at(i) = scaL.at(i);
            matMult(-1, scaA.at(i), dx.row(normals.blockIndex(i), normals.blockSize(i)), e.at(i));

            // redundancy
            A.at(i) = scaA.at(i);
            redundancy.at(i) = computeRedundancy(A, normals);
            A.at(i) = Matrix();
          } // for(i=epochCount)

        estimateAccuracy(e, redundancy, perAxisSca, huber, huberPower, sigmaAxisSca, sigmaEpochSca);
      }

      // Accuracy & outlier detection (accelerometer)
      // --------------------------------------------
      {
        std::vector<Vector> e(epochCount);
        std::vector<Vector> redundancy(epochCount);
        std::vector<Matrix> A(epochCount+1);
        for(UInt i=0; i<epochCount; i++)
          if(accL.at(i).size())
          {
            // decorrelated residuals
            e.at(i) = accL.at(i);
            for(UInt k=0; k<accA.at(i).size(); k++)
              matMult(-1, accA.at(i).at(k), dx.row(normals.blockIndex(accIdxA.at(i)+k), normals.blockSize(accIdxA.at(i)+k)), e.at(i));
            if(countAccCal)
              matMult(-1., accB.at(i), dx.row(normals.blockIndex(idxBlockCal), normals.blockSize(idxBlockCal)), e.at(i));

            // redundancy
            for(UInt k=0; k<accA.at(i).size(); k++)
              A.at(accIdxA.at(i)+k) = accA.at(i).at(k);
            A.back() = accB.at(i);
            redundancy.at(i) = computeRedundancy(A, normals);
            for(UInt k=0; k<accA.at(i).size(); k++)
              A.at(accIdxA.at(i)+k) = Matrix();
          } // for(i=epochCount)

        estimateAccuracy(e, redundancy, perAxisAcc, huber, huberPower, sigmaAxisAcc, sigmaEpochAcc);
      }

      // convergence criteria
      // --------------------
      if(fabs(sigma0-sigma0Old)<1e-5)
        break;
    } // for(iter)

    //=================================

    // full variance covariance matrix
    Matrix C;
    if(!fileNameOutCovariance.empty() || !fileNameOutCovarianceMatrix.empty())
    {
      C = Matrix(normals.blockIndex(blockCount), Matrix::TRIANGULAR, Matrix::UPPER);
      for(UInt i=0; i<blockCount; i++)
        for(UInt k=i; k<blockCount; k++)
          if(normals.isBlockUsed(i, k))
            copy(normals.N(i,k), C.slice(normals.blockIndex(i), normals.blockIndex(k), normals.blockSize(i), normals.blockSize(k)));
      cholesky2Inverse(C);
      C *= sigma0 * sigma0;
    }

    if(!fileNameOutCovarianceMatrix.empty())
      writeFileMatrix(fileNameOutCovarianceMatrix.appendBaseName(".arc"+arcNo%"%03i"s), C.slice(0,0,3*epochCount,3*epochCount));

    // Results
    ArcObservation arc;
    arc.sigmaSca = sigmaAxisSca;
    arc.sigmaAcc = sigmaAxisAcc;
    arc.x        = xAccCal;

    // star camera arc
    if(!fileNameOutStarCamera.empty())
      for(UInt i=0; i<epochCount; i++)
        if(hasAcc.at(i) || hasSca.at(i))
        {
          StarCameraEpoch epoch;
          epoch.time   = times.at(i);
          epoch.rotary = Rotary3d(q.row(i).trans());
          arc.starCamera.push_back(epoch);
        }

    // star camera covariance arc
    if(!fileNameOutCovariance.empty())
    {
      for(UInt i=0; i<epochCount; i++)
        if(hasAcc.at(i) || hasSca.at(i))
        {
          Covariance3dEpoch epoch;
          epoch.time       = times.at(i);
          epoch.covariance = Tensor3d(C.slice(3*i,3*i,3,3));
          arc.starCameraCovariance.push_back(epoch);
        }
    }

    // angular acceleration observation (calibrated)
    if(!fileNameOutAngularAcc.empty())
    {
      UInt idxAcc = 0;
      for(UInt i=0; i<epochCount; i++)
        if(index(times.at(i), timesAcc, idxAcc))
          arc.angularAcc.push_back(angularAcc.at(idxAcc));
    }

    // star Camera epoch sigmas
    if(!fileNameOutEpochSigmasStarCamera.empty())
    {
      for(UInt i=0; i<epochCount; i++)
        if(hasSca.at(i))
        {
          MiscValueEpoch epoch;
          epoch.time  = times.at(i);
          epoch.value = sigmaEpochSca(i);
          arc.epochSigmaSca.push_back(epoch);
        }
    }

    // accelerometer epoch sigmas
    if(!fileNameOutEpochSigmasAccelerometer.empty())
    {
      for(UInt i=0; i<epochCount; i++)
        if(hasAcc.at(i))
        {
          MiscValueEpoch epoch;
          epoch.time  = times.at(i);
          epoch.value = sigmaEpochAcc(i);
          arc.epochSigmaAcc.push_back(epoch);
        }
    }

    return arc;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector InstrumentStarCameraAngularAccelerometerFusion::computeRedundancy(std::vector<Matrix> &A, MatrixDistributed &normals, Double threshold)
{
  try
  {
    Vector redundancy = {1., 1., 1.};
    for(UInt i=0; i<normals.blockCount()-1; i++)
      if(A.at(i).size())
      {
        triangularSolve(1., normals.N(i,i).trans(), A.at(i).trans());
        for(UInt idAxis=0; idAxis<3; idAxis++)
          redundancy(idAxis) -= quadsum(A.at(i).row(idAxis));
        for(UInt k=i+1; k<normals.blockCount(); k++)
          if(normals.isBlockUsed(i,k))
          {
            if(!A.at(k).size())
              A.at(k) = Matrix(A.at(i).rows(), normals.blockSize(k));
            matMult(-1., A.at(i), normals.N(i,k), A.at(k));
          }
        if(quadsum(A.at(i)) < threshold)
          break;
      }
    if(normals.isBlockUsed(normals.blockCount()-1, normals.blockCount()-1) && A.back().size())
    {
      triangularSolve(1., normals.N(normals.blockCount()-1,normals.blockCount()-1).trans(), A.back().trans());
      for(UInt idAxis=0; idAxis<3; idAxis++)
        redundancy(idAxis) -= quadsum(A.back().row(idAxis));
    }

    return redundancy;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void InstrumentStarCameraAngularAccelerometerFusion::estimateAccuracy(std::vector<Vector> &e, std::vector<Vector> &redundancy,
                                                                      Bool sigmaPerAxis, Double huber, Double huberPower,
                                                                      Vector &sigmaAxis, Vector &sigmaEpoch)
{
  try
  {
    // total accuracy
    Double ePe_sum       = 0;
    Double redundany_sum = 0;
    for(UInt idAxis=0; idAxis<3; idAxis++)
    {
      for(UInt i=0; i<e.size(); i++)
        if(redundancy.at(i).size() && (sigmaEpoch(i) <= 1))
        {
          ePe_sum       += std::pow(e.at(i)(idAxis)/sigmaEpoch(i), 2);
          redundany_sum += redundancy.at(i)(idAxis);
        }

      if(sigmaPerAxis)
      {
        Double factor = Vce::standardDeviation(ePe_sum, redundany_sum, huber, huberPower);
        for(UInt i=0; i<e.size(); i++)
          if(redundancy.at(i).size())
            e.at(i)(idAxis) *= 1./factor;
        sigmaAxis(idAxis) *= factor;
        ePe_sum = redundany_sum = 0;
      }
    }

    if(!sigmaPerAxis)
    {
      Double factor = Vce::standardDeviation(ePe_sum, redundany_sum, huber, huberPower);
      for(UInt i=0; i<e.size(); i++)
        if(redundancy.at(i).size())
          e.at(i) *= 1/factor;
      for(UInt idAxis=0; idAxis<3; idAxis++)
        sigmaAxis(idAxis) *= factor;
    }

    // Outlier detection
    for(UInt i=0; i<e.size(); i++)
      if(redundancy.at(i).size())
      {
        const Double s = sigmaEpoch(i)*std::sqrt(quadsum(e.at(i))/(3*mean(redundancy.at(i))));
        sigmaEpoch(i) = 1;
        if((s > huber) && (mean(redundancy.at(i))>1e-3))
          sigmaEpoch(i) = std::pow(s/huber, huberPower);
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
