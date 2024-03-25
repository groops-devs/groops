/***********************************************/
/**
* @file instrumentAccelerometerEstimateParameters.cpp
*
* @brief Estimate parameters from accelerometer data.
*
* @author Andreas Kvas
* @date 2022-07-27
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program estimates calibration parameters for acceleration data given given an optional reference acceleration.
Specifically, the program solves the equation
\begin{equation}
  \mathbf{a} - \mathbf{a}_\text{ref} = \mathbf{f}(\mathbf{x}) + \mathbf{e}
\end{equation}
for the unknown parameters $\mathbf{x}$, where $\mathbf{a}$ is given in \configFile{inputfileAccelerometer}{instrument} and
$\mathbf{a}_\text{ref}$ is given in  \configFile{inputfileAccelerometerReference}{instrument}.
The parametrization of $\mathbf{x}$ can be set via \configClass{parametrizationAcceleration}{parametrizationAccelerationType}.
Optionally, the empirical covariance functions for the accelerations $\mathbf{a}$ can be estimated by enabling \config{estimateCovarianceFunctions}.

The estimated parameters are written to the file \configFile{outputfileSolution}{matrix} and can be used by
\program{InstrumentAccelerometerApplyEstimatedParameters} to calibrate accelerometer measurements.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"
#include "files/fileSatelliteModel.h"
#include "files/fileParameterName.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/ephemerides/ephemerides.h"
#include "classes/parametrizationAcceleration/parametrizationAcceleration.h"
#include "misc/varianceComponentEstimation.h"

/***** CLASS ***********************************/

/** @brief  Estimate parameters from accelerometer data.
* @ingroup programsGroup */
class InstrumentAccelerometerEstimateParameters
{
  InstrumentFile                 accFile, accFileSim, orbitFile, starCameraFile;
  EarthRotationPtr               earthRotation;
  EphemeridesPtr                 ephemerides;
  ParametrizationAccelerationPtr parameterAcceleration;
  SatelliteModelPtr              satellite;

  Matrix N;        // normal equation matrix
  Vector n;        // right hand sides
  Vector x;        // solution
  Matrix Wz;       // monte carlo vector for redundancy computation
  Double lPl;      // =l'Pl, weighted norm of the observations
  UInt   obsCount; // number of observations

  Bool   estimateCovarianceFunctionVCE;
  Bool   estimateArcSigmas;
  Bool   estimateEpochSigmas;

  Double sampling;
  Matrix covFunc, CosTransform, Psd;
  Vector arcSigmas, arcSigmasNew;
  Matrix ePe, redundancy;  // one row for each frequency, one column for each component
  std::vector<ObservationSigmaArc> arcListEpochSigma;

  class ObservationEquationArc
  {
    public:
      std::vector<Time> times;
      Vector l;
      Matrix A, B;
      Vector xArc; // arc wise parameters
  };

  std::vector<ObservationEquationArc> observationEquations;

  void computeObservationEquation(UInt arcNo);
  void decorrelate(const std::vector<Time> &times, Double sigmaArc, const ObservationSigmaArc &sigmaEpoch, std::vector<Matrix> &W, const std::list<MatrixSlice> &A);
  void buildNormals(UInt arcNo);
  void computeRedundancies(UInt arcNo);
  ObservationSigmaArc computeEpochSigmas(UInt arcNo);

public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(InstrumentAccelerometerEstimateParameters, PARALLEL, "estimate accelerometer parameters.", Instrument)

/***********************************************/

void InstrumentAccelerometerEstimateParameters::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    estimateCovarianceFunctionVCE = estimateArcSigmas = estimateEpochSigmas = FALSE;

    FileName fileNameOutSolution, fileNameOutParameterName;
    FileName fileNameOutCovFunc, fileNameOutArcSigmas, fileNameOutEpochSigmas;
    FileName fileNameInAccelerometer,  fileNameInAccelerometerReference, fileNameInOrbit, fileNameInStarCamera, fileNameInSatelliteModel;
    Double   sigmaX, sigmaY, sigmaZ;
    UInt     maxIter;

    readConfig(config, "outputfileSolution",              fileNameOutSolution,        Config::MUSTSET,  "",     "values for estimated parameters");
    readConfig(config, "outputfileParameterNames",        fileNameOutParameterName,   Config::OPTIONAL, "",     "names of the estimated parameters");
    if(readConfigSequence(config, "estimateArcSigmas", Config::OPTIONAL, "", ""))
    {
      estimateArcSigmas = TRUE;
      readConfig(config, "outputfileArcSigmas",           fileNameOutArcSigmas,   Config::OPTIONAL, "", "accuracies of each arc");
      endSequence(config);
    }
    if(readConfigSequence(config, "estimateEpochSigmas", Config::OPTIONAL, "", ""))
    {
      estimateEpochSigmas = TRUE;
      readConfig(config, "outputfileEpochSigmas",         fileNameOutEpochSigmas, Config::OPTIONAL, "", "estimated epoch-wise sigmas");
      endSequence(config);
    }
    if(readConfigSequence(config, "estimateCovarianceFunctions", Config::OPTIONAL, "", ""))
    {
      estimateCovarianceFunctionVCE = TRUE;
      readConfig(config, "outputfileCovarianceFunction",  fileNameOutCovFunc,     Config::OPTIONAL, "", "covariance functions for x, y, z direction");
      endSequence(config);
    }
    readConfig(config, "inputfileAccelerometer",          fileNameInAccelerometer,          Config::MUSTSET,  "",      "");
    readConfig(config, "inputfileAccelerometerReference", fileNameInAccelerometerReference, Config::OPTIONAL, "",      "if not given, reference acceleration is assumed zero");
    readConfig(config, "inputfileOrbit",                  fileNameInOrbit,                  Config::OPTIONAL, "",      "may be needed by parametrizationAcceleration");
    readConfig(config, "inputfileStarCamera",             fileNameInStarCamera,             Config::OPTIONAL, "",      "may be needed by parametrizationAcceleration");
    readConfig(config, "inputfileSatelliteModel",         fileNameInSatelliteModel,         Config::OPTIONAL, "{groopsDataDir}/satelliteModel/", "satellite macro model (may be needed by parametrizationAcceleration)");
    readConfig(config, "earthRotation",                   earthRotation,                    Config::OPTIONAL, "file",  "may be needed by parametrizationAcceleration");
    readConfig(config, "ephemerides",                     ephemerides,                      Config::OPTIONAL, "jpl",   "may be needed by parametrizationAcceleration");
    readConfig(config, "parametrizationAcceleration",     parameterAcceleration,            Config::MUSTSET,  "",      "");
    readConfig(config, "sigmaX",                          sigmaX,                           Config::DEFAULT,  "1e-9",  "apriori accuracy in x-axis");
    readConfig(config, "sigmaY",                          sigmaY,                           Config::DEFAULT,  "10e-9", "apriori accuracy in y-axis");
    readConfig(config, "sigmaZ",                          sigmaZ,                           Config::DEFAULT,  "1e-9",  "apriori accuracy in z-axis");
    readConfig(config, "iterationCount",                  maxIter,                          Config::DEFAULT,  "5",     "iteration count for determining the covariance function");
    if(isCreateSchema(config)) return;

     // ======================================================

    // init instrument files + satellite model
    // ---------------------------------------
    logStatus<<"read instrument data <"<<fileNameInAccelerometer<<">"<<Log::endl;
    logStatus<<"read reference data <"<<fileNameInAccelerometerReference<<">"<<Log::endl;
    accFile.open(fileNameInAccelerometer);
    accFileSim.open(fileNameInAccelerometerReference);
    orbitFile.open(fileNameInOrbit);
    starCameraFile.open(fileNameInStarCamera);
    InstrumentFile::checkArcCount({accFile, accFileSim, orbitFile, starCameraFile});
    const UInt arcCount = accFile.arcCount();

    if(!fileNameInSatelliteModel.empty())
    {
      logStatus<<"read satellite model <"<<fileNameInSatelliteModel<<">"<<Log::endl;
      readFileSatelliteModel(fileNameInSatelliteModel, satellite);
    }

    if(Parallel::isMaster(comm) && !fileNameOutParameterName.empty())
    {
      logStatus<<"write parameter names to <"<<fileNameOutParameterName<<">"<<Log::endl;
      std::vector<ParameterName> names, nameArc;
      parameterAcceleration->parameterName(names);
      parameterAcceleration->parameterNameArc(nameArc);
      for(UInt arcNo=0; arcNo<arcCount; arcNo++)
      {
        const std::string str = "arc"+arcNo%"%i"s+".";
        for(UInt i=0; i<nameArc.size(); i++)
        {
          ParameterName param = nameArc.at(i);
          param.type = str+param.type;
          names.push_back(param);
        }
      }
      writeFileParameterName(fileNameOutParameterName, names);
    }

    // ======================================================

    // setup observation equations
    // ---------------------------
    logStatus<<"compute observation equations"<<Log::endl;
    observationEquations.resize(arcCount);
    const std::vector<UInt> processNo = Parallel::forEach(arcCount, [&](UInt arcNo) {computeObservationEquation(arcNo);}, comm);

    // determine sampling
    sampling  = 0.;
    for(const auto &arc : observationEquations)
      if(arc.times.size())
        sampling += medianSampling(arc.times).seconds();
    Parallel::reduceSum(sampling, 0, comm);
    Parallel::broadCast(sampling, 0, comm);
    sampling /= arcCount;

    // Determine max. length of ovariance functions
    UInt covLength = 0;
    for(const auto &arc : observationEquations)
      if(arc.times.size())
        covLength = std::max(covLength, static_cast<UInt>(std::round((arc.times.back()-arc.times.front()).seconds()/sampling)+1));
    Parallel::reduceMax(covLength, 0, comm);
    Parallel::broadCast(covLength, 0, comm);

    logInfo<<"  length of covariance function: "<<covLength<<" epochs with a sampling of "<<sampling<<" seconds"<<Log::endl;

    // init arc sigmas and covariance function
    // ---------------------------------------
    arcSigmas    = Vector(arcCount, 1.0);
    CosTransform = Vce::cosTransform(covLength);
    covFunc = Vce::readCovarianceFunction(FileName(), covLength, 3, sampling);
    covFunc.column(1)  *= std::pow(sigmaX, 2);
    covFunc.column(2)  *= std::pow(sigmaY, 2);
    covFunc.column(3)  *= std::pow(sigmaZ, 2);
    Psd = CosTransform * covFunc.column(1, 3);

    // init epoch sigmas
    // -----------------
    arcListEpochSigma.resize(arcCount);
    Parallel::forEachProcess(arcCount, [this](UInt arcNo)
    {
      const auto &times = observationEquations.at(arcNo).times;
      arcListEpochSigma.at(arcNo) = Arc(times, Matrix(times.size(), 2), Epoch::OBSERVATIONSIGMA);
    }, processNo, comm, FALSE/*timing*/);

    // start iteration
    // ---------------
    const UInt countParameter = parameterAcceleration->parameterCount();
    for(UInt iter=0; iter<maxIter; iter++)
    {
      if(maxIter>1) logStatus<<"starting iteration "<<iter+1<<Log::endl;

      // solve normal equations
      // ----------------------
      if(countParameter)
      {
        logStatus<<"accumulate system of normal equations"<<Log::endl;
        N        = Matrix(countParameter, Matrix::SYMMETRIC);
        n        = Vector(countParameter);
        lPl      = 0;
        obsCount = 0;
        Parallel::forEachProcess(arcCount, [this](UInt arcNo) {buildNormals(arcNo);}, processNo, comm);
        Parallel::reduceSum(N,        0, comm);
        Parallel::reduceSum(n,        0, comm);
        Parallel::reduceSum(obsCount, 0, comm);
        Parallel::reduceSum(lPl,      0, comm);

        logStatus<<"solve system of normal equations"<<Log::endl;
        if(Parallel::isMaster(comm))
        {
          // regularize unused parameters
          UInt countUnused = 0;
          for(UInt k=0; k<N.rows(); k++)
            if(N(k, k) == 0.0)
            {
              N(k, k) = 1.0;
              countUnused++;
            }
          obsCount += countUnused;
          if(countUnused)
            logInfo<<"  "<<countUnused<<" parameters unused -> set to zero"<<Log::endl;

          x = solve(N, n);
          logInfo<<"  aposteriori sigma = "<<std::sqrt((lPl-inner(x, n))/(obsCount-x.rows()))<<Log::endl;

          Wz = Vce::monteCarlo(x.rows(), 100); // monte carlo vector for VCE
          triangularSolve(1., N, Wz);
        }
        Parallel::broadCast(x,  0, comm);
        Parallel::broadCast(Wz, 0, comm);
      }
      Parallel::barrier(comm);

      // compute redundancies
      // --------------------
      logStatus<<"compute arc parameters and redundancies"<<Log::endl;
      arcSigmasNew = Vector(arcCount);
      ePe = redundancy = Matrix(covLength, 3);
      Parallel::forEachProcess(arcCount, [this](UInt arcNo) {computeRedundancies(arcNo);}, processNo, comm);
      Parallel::barrier(comm);

      // write parameter vector
      // ----------------------
      if(!fileNameOutSolution.empty())
      {
        // collect solution vector
        std::vector<Vector> solutions(arcCount);
        Parallel::forEachProcess(solutions, [this](UInt arcNo) {return observationEquations.at(arcNo).xArc;}, processNo, comm, FALSE/*timing*/);

        if(Parallel::isMaster(comm))
        {
          // copy global and arc parameters to one vector
          solutions.insert(solutions.begin(), x); // global parametes
          const UInt count = std::accumulate(solutions.begin(), solutions.end(), UInt(0), [](UInt c, const Vector &x) {return c+x.size();});

          Vector xOut(count);
          UInt   idx = 0;
          for(const Vector &solution : solutions)
          {
            copy(solution, xOut.row(idx, solution.rows()));
            idx += solution.rows();
          }

          logStatus<<"write solution to <"<<fileNameOutSolution<<">"<<Log::endl;
          writeFileMatrix(fileNameOutSolution, xOut);
        }
      }

      // sigmas per arc
      // --------------
      if(estimateArcSigmas)
      {
        Parallel::reduceSum(arcSigmasNew, 0, comm);
        if(Parallel::isMaster(comm))
        {
          arcSigmas = arcSigmasNew;
          const Double sigma0 = Vce::meanSigma(arcSigmas);
          arcSigmas *= 1./sigma0;
          logInfo<<"  sigma per arc (median): "<<sigma0<<Log::endl;

          if(!fileNameOutArcSigmas.empty())
          {
            logStatus<<"write arc sigma file <"<<fileNameOutArcSigmas<<">"<<Log::endl;
            writeFileMatrix(fileNameOutArcSigmas, arcSigmas);
          }
        }
        Parallel::broadCast(arcSigmas, 0, comm);
      }

      // sigmas per epoch
      // --------------
      if(estimateEpochSigmas)
      {
        Parallel::forEachProcess(arcListEpochSigma, [this](UInt arcNo) {return computeEpochSigmas(arcNo);}, processNo, comm, FALSE/*timing*/);
        if(Parallel::isMaster(comm))
        {
          UInt count = std::accumulate(arcListEpochSigma.begin(), arcListEpochSigma.end(), UInt(0), [](UInt c, auto &arc) {return c+arc.size();});
          UInt countOutlier = 0;
          for(auto &arc : arcListEpochSigma)
            countOutlier += std::count_if(arc.begin(), arc.end(), [](const ObservationSigmaEpoch &e) {return (e.sigma > 0);});
          logInfo<<"  "<<countOutlier<<" of "<<count<<" outliers ("<<100.*countOutlier/count%"%.2f%%)"s<<Log::endl;
        }
        if(Parallel::isMaster(comm) && !fileNameOutEpochSigmas.empty())
        {
          logStatus<<"write epoch sigma file <"<<fileNameOutEpochSigmas<<">"<<Log::endl;
          InstrumentFile::write(fileNameOutEpochSigmas, arcListEpochSigma);
        }
      }

      // estimate new PSD through variance component estimation
      // ------------------------------------------------------
      if(estimateCovarianceFunctionVCE)
      {
        logStatus<<"compute psd through variance component estimation"<<Log::endl;
        Parallel::reduceSum(ePe, 0, comm);
        Parallel::reduceSum(redundancy, 0, comm);
        if(Parallel::isMaster(comm))
        {
          Double maxFactor = 0;
          Vce::estimatePsd(ePe, redundancy, Psd, maxFactor);
          logInfo<<"  max. PSD adjustment factor: "<<maxFactor<<Log::endl;
        }
        Parallel::broadCast(Psd, 0, comm);
        copy(CosTransform*Psd, covFunc.column(1, Psd.columns())); // compute new covariance function

        if(Parallel::isMaster(comm) && !fileNameOutCovFunc.empty())
        {
          logStatus<<"write covariance function file <"<<fileNameOutCovFunc<<">"<<Log::endl;
          writeFileMatrix(fileNameOutCovFunc, covFunc);
        }
      }
    } // for(iter)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void InstrumentAccelerometerEstimateParameters::computeObservationEquation(UInt arcNo)
{
  try
  {
    AccelerometerArc acc        = accFile.readArc(arcNo);
    AccelerometerArc accSim     = accFileSim.readArc(arcNo);
    OrbitArc         orbit      = orbitFile.readArc(arcNo);
    StarCameraArc    starCamera = starCameraFile.readArc(arcNo);
    Arc::checkSynchronized({acc, accSim, orbit, starCamera});

    const std::vector<Time> times = acc.times();
    parameterAcceleration->setIntervalArc(times.front(), times.back()+medianSampling(times));

    const UInt epochCount = acc.size();
    Vector l(3*epochCount);
    Matrix A(3*epochCount, parameterAcceleration->parameterCount());
    Matrix B(3*epochCount, parameterAcceleration->parameterCountArc());

    for(UInt i=0; i<epochCount; i++)
    {
      Rotary3d rotEarth, rotSat;
      Vector3d position, velocity;
      if(earthRotation)     rotEarth = earthRotation->rotaryMatrix(acc.at(i).time);
      if(starCamera.size()) rotSat   = starCamera.at(i).rotary;
      if(orbit.size())      position = orbit.at(i).position;
      if(orbit.size())      velocity = orbit.at(i).velocity;

      Vector lk = acc.at(i).acceleration.vector();
      if(accSim.size())
        lk -= accSim.at(i).acceleration.vector();
      Matrix Ak(3, parameterAcceleration->parameterCount());
      Matrix Bk(3, parameterAcceleration->parameterCountArc());
      parameterAcceleration->compute(satellite, acc.at(i).time, position, velocity, rotSat, rotEarth, ephemerides, Ak, Bk);

      const Matrix R = (rotEarth*rotSat).matrix();  // rotate into accelerometer frame
      l(i)              = lk(0);
      l(i+epochCount)   = lk(1);
      l(i+2*epochCount) = lk(2);
      if(A.size())
      {
        Ak = R.trans() * Ak;
        copy(Ak.row(0), A.row(i));
        copy(Ak.row(1), A.row(i+epochCount));
        copy(Ak.row(2), A.row(i+2*epochCount));
      }
      if(B.size())
      {
        Bk = R.trans() * Bk;
        copy(Bk.row(0), B.row(i));
        copy(Bk.row(1), B.row(i+epochCount));
        copy(Bk.row(2), B.row(i+2*epochCount));
      }
    }

    observationEquations.at(arcNo).l = l;
    observationEquations.at(arcNo).A = A;
    observationEquations.at(arcNo).B = B;
    observationEquations.at(arcNo).times = times;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void InstrumentAccelerometerEstimateParameters::decorrelate(const std::vector<Time> &times, Double sigmaArc, const ObservationSigmaArc &sigmaEpoch,
                                                            std::vector<Matrix> &W, const std::list<MatrixSlice> &A)
{
  try
  {
    std::vector<UInt> index;
    for(const Time &t : times)
      index.push_back(static_cast<UInt>(std::round((t-times.front()).seconds()/sampling)));

    W = std::vector<Matrix>(3, Matrix(times.size(), Matrix::SYMMETRIC, Matrix::UPPER));
    for(UInt idAxis=0; idAxis<3; idAxis++)
    {
      for(UInt i=0; i<times.size(); i++)
        for(UInt k=i; k<times.size(); k++)
          W.at(idAxis)(i, k) = sigmaArc*sigmaArc * covFunc(index.at(k-i), idAxis+1);

      if(sigmaEpoch.size())
        for(UInt i=0; i<times.size(); i++)
            W.at(idAxis)(i, i) += std::pow(sigmaEpoch.at(i).sigma, 2);

      cholesky(W.at(idAxis));
      for(MatrixSliceRef WA : A)
        if(WA.size())
          triangularSolve(1., W.at(idAxis).trans(), WA.row(idAxis*times.size(), times.size()));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void InstrumentAccelerometerEstimateParameters::buildNormals(UInt arcNo)
{
  try
  {
    // decorrelate
    Vector Wl = observationEquations.at(arcNo).l;
    Matrix WA = observationEquations.at(arcNo).A;
    Matrix WB = observationEquations.at(arcNo).B;
    std::vector<Matrix> W;
    decorrelate(observationEquations.at(arcNo).times, arcSigmas(arcNo), arcListEpochSigma.at(arcNo), W, {Wl, WA, WB});

    // estimate arc dependent parameters
    Vector tau;
    if(WB.size())
    {
      tau = QR_decomposition(WB);
      QTransMult(WB, tau, Wl); // transform observations: l:= Q'l
      QTransMult(WB, tau, WA); // transform design matrix A:=Q'A
    }
    // use only nullspace of design matrix WB
    MatrixSlice A_bar( WA.row(WB.columns(), WA.rows()-WB.columns()) );
    MatrixSlice l_bar( Wl.row(WB.columns(), Wl.rows()-WB.columns()) );

    // build normals
    this->lPl      += quadsum(l_bar);
    this->obsCount += l_bar.rows();
    matMult(1., A_bar.trans(), l_bar, n);
    rankKUpdate(1., A_bar, N);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void InstrumentAccelerometerEstimateParameters::computeRedundancies(UInt arcNo)
{
  try
  {
    Vector We = observationEquations.at(arcNo).l;
    Matrix WA = observationEquations.at(arcNo).A;
    Matrix WB = observationEquations.at(arcNo).B;
    if(WA.size())
      matMult(-1., WA,  x, We);

    std::vector<Matrix> W;
    decorrelate(observationEquations.at(arcNo).times, arcSigmas(arcNo), arcListEpochSigma.at(arcNo), W, {We, WA, WB});

    // estimate & eliminate arc dependent parameters
    // ---------------------------------------------
    if(WB.size())
    {
      Vector tau = QR_decomposition(WB);
      QTransMult(WB, tau, We);           // transform observations: l:= Q'l
      observationEquations.at(arcNo).xArc = We.row(0, tau.rows());
      triangularSolve(1., WB.row(0, tau.rows()), observationEquations.at(arcNo).xArc);
      We.row(0, WB.columns()).setNull(); // residuals: remove WB*x
      QMult(WB, tau, We);                // back transformation

      if(WA.size())
      {
        QTransMult(WB, tau, WA);           // transform design matrix A:=Q'A
        WA.row(0, WB.columns()).setNull(); // residuals: remove WB*x
        QMult(WB, tau, WA);                // back transformation
      }

      generateQ(WB, tau);
      WB.setType(Matrix::GENERAL);
    }

    if(!estimateArcSigmas && !estimateCovarianceFunctionVCE)
      return;

    // Without variance component estimation
    // -------------------------------------
    Matrix WAz(We.rows(), Wz.columns());
    if(WA.size())
      matMult( 1., WA, Wz, WAz);
    if(!estimateCovarianceFunctionVCE)
    {
      const Double redundancy = We.rows() - quadsum(WAz) - quadsum(WB);
      arcSigmasNew(arcNo) = std::sqrt(quadsum(We)/redundancy) * arcSigmas(arcNo);
      return;
    }

    // Variance component estimation
    // -----------------------------
    std::vector<UInt> index;
    for(const Time &t : observationEquations.at(arcNo).times)
      index.push_back(static_cast<UInt>(std::round((t-observationEquations.at(arcNo).times.front()).seconds()/sampling)));
    const UInt countEpoch = observationEquations.at(arcNo).times.size();
    Double ePeSum=0, redundancySum=0;

    for(UInt idAxis=0; idAxis<Psd.columns(); idAxis++)
    {
      Matrix R;
      Vector WWe;
      Vce::redundancy(W.at(idAxis), We.row(idAxis*countEpoch, countEpoch), WAz.row(idAxis*countEpoch, countEpoch), WB.row(idAxis*countEpoch, countEpoch), R, WWe);
      Vce::psd(R, WWe, index, arcSigmas(arcNo), CosTransform, Psd.column(idAxis),
               ePe.column(idAxis), redundancy.column(idAxis), ePeSum, redundancySum);
    }
    arcSigmasNew(arcNo) = std::sqrt(ePeSum/redundancySum) * arcSigmas(arcNo);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

ObservationSigmaArc InstrumentAccelerometerEstimateParameters::computeEpochSigmas(UInt arcNo)
{
  try
  {
    const Double huber = 2.5;
    const Double threshold2 = std::pow(huber * arcSigmas(arcNo), 2) * sum(covFunc.slice(0, 1, 1, 3));
    const UInt   epochCount = observationEquations.at(arcNo).times.size();

    Vector e = observationEquations.at(arcNo).l;
    if(observationEquations.at(arcNo).A.size())
      matMult(-1.0, observationEquations.at(arcNo).A, x, e);
    if(observationEquations.at(arcNo).B.size())
      matMult(-1, observationEquations.at(arcNo).B, observationEquations.at(arcNo).xArc, e);

    for(UInt i=0; i<epochCount; i++)
    {
      const Double e2 = quadsum(Vector({e(i), e(epochCount+i), e(2*epochCount+i)}));
      arcListEpochSigma.at(arcNo).at(i).sigma = 0.;
      if(e2 > threshold2)
        arcListEpochSigma.at(arcNo).at(i).sigma = std::sqrt((e2-threshold2)/3);
    }

    return arcListEpochSigma.at(arcNo);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
