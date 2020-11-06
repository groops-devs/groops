/***********************************************/
/**
* @file preprocessingPod.cpp
*
* @brief Estimate covariance function / arc weights.
*
* @author Torsten Mayer-Guerr
* @date 2011-06-30
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program estimates empirical covariance functions of the instrument noise and determine arc wise variances to
downweight arc with outliers.

A complete least squares adjustment for gravity field determination is performed by computing the \config{observation}
equations, see \configClass{observation:podIntegral}{observationType:podIntegral} or
\configClass{observation:podVariational}{observationType:podVariational} for details. The normal equations
are accumulated and solved to \configFile{outputfileSolution}{matrix} together with the estimated accuracies
\configFile{outputfileSigmax}{matrix}. The estimated residuals~$\hat{\M e}=\M l-\M A\hat{\M x}$ can be computed with
\config{computeResiduals}.

For each component (along, cross, radial) of the kinematic orbit positions a noise covariance function is estimated
\begin{equation}
  \text{cov}(\Delta t_i) = \sum_{n=0}^{N-1} a_n^2 \cos\left(\frac{\pi}{T} n\Delta t_i\right).
\end{equation}
The covariance matrix is composed of the sum of matrices $F_n$ and unknown variance factors
\begin{equation}
  \M\Sigma = a_1^2\M F_1 + a_2^2 \M F_2 + \cdots + a_N^2\M F_N,
\end{equation}
with the cosine transformation matrices
\begin{equation}
  \M F_n = \left(\cos\left(\frac{\pi}{T} n(t_i-t_k)\right)\right)_{ik}.
\end{equation}

An additional variance factor can be computed (\config{estimateArcSigmas}) for each arc~$k$  according to
\begin{equation}
  \hat{\sigma}_k^2 = \frac{\hat{\M e}_k^T\M\Sigma^{-1}\hat{\M e}_k}{r_k},
\end{equation}
where $r_k$ is the redundancy. This variance factor should be around one for normal behaving arcs
as the noise characteristics is already considered by the covariance matrix but bad arcs get a much larger variance.
By appling this factor bad arcs or arcs with large outliers are downweighted.
)";

/***********************************************/

#include "programs/program.h"
#include "parallel/matrixDistributed.h"
#include "files/fileArcList.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"
#include "files/fileParameterName.h"
#include "misc/observation/observationMiscPod.h"
#include "misc/observation/covariancePod.h"
#include "misc/varianceComponentEstimation.h"

/***** CLASS ***********************************/

/** @brief Estimate covariance function / arc weights.
* @ingroup programsGroup */
class PreprocessingPod
{
public:
  void run(Config &config);

private:
  ObservationMiscPodPtr observationMisc;
  std::vector<ObservationMiscPod::Arc> observationArc;
  InstrumentFile fileCovPodEpoch;

  std::vector< std::vector<UInt> > arcIndexN;
  std::vector< std::vector<UInt> > arcIndexA;

  Bool estimateCovarianceFunctionVCE;
  Bool estimateArcSigmas;
  Bool estimateResiduals;

  // normal equations
  // ----------------
  MatrixDistributed normals;
  Matrix n;        // right hand sides
  Matrix x;        // solution
  Matrix Wz;       // monte carlo vector for redundancy computation
  Double lPl;      // =l'Pl, weighted norm of the observations
  UInt   obsCount; // number of observations

  // covariance
  // ----------
  Vector  sigmaPod, sigmaPodNew;
  Double  samplingPod;
  Matrix  covFuncPod;         // (covLength x 3) one column for x,y,z
  Matrix  PsdPod;
  Matrix  ePePod;
  Matrix  redundancyPod;      // one row for each frequency, one column for each component
  Matrix  CosTransformPod;

  // residuals
  std::vector<Arc> arcListPod;

  void computeObservationEquation(UInt arcNo);
  void buildNormals(UInt arcNo);
  void computeRedundancies(UInt arcNo);
  void computeResiduals(UInt arcNo);
  Arc  collectArcPod(UInt arcNo) const {return arcListPod.at(arcNo);}
};

GROOPS_REGISTER_PROGRAM(PreprocessingPod, PARALLEL, "Estimate covariance function / arc weights.", Preprocessing, Covariance, Residuals)

/***********************************************/

void PreprocessingPod::run(Config &config)
{
  try
  {
    estimateCovarianceFunctionVCE = estimateArcSigmas = estimateResiduals = FALSE;

    FileName fileNameSolution, fileNameSigmax, fileNameParaName;
    FileName fileNameArcSigmaPod, fileNameCovPod, fileNameResidualsPod;
    FileName fileNameArcList;
    FileName sigmaPodName, covPodEpochName, covPodName;
    Double   sigma0Pod=1;
    Double   adjustmentThreshold;
    UInt     iterCount;

    renameDeprecatedConfig(config, "arcList", "inputfileArcList", date2time(2020, 7, 7));

    readConfig(config, "outputfileSolution",      fileNameSolution, Config::OPTIONAL, "", "estimated parameter vector (static part only)");
    readConfig(config, "outputfileSigmax",        fileNameSigmax,   Config::OPTIONAL, "", "standard deviations of the parameters (sqrt of the diagonal of the inverse normal equation)");
    readConfig(config, "outputfileParameterName", fileNameParaName, Config::OPTIONAL, "", "names of estimated parameters (static part only)");
    if(readConfigSequence(config, "estimateArcSigmas", Config::OPTIONAL, "", ""))
    {
      estimateArcSigmas = TRUE;
      readConfig(config, "outputfileSigmasPerArcPod", fileNameArcSigmaPod, Config::OPTIONAL, "", "accuracies of each arc (POD2)");
      endSequence(config);
    }
    if(readConfigSequence(config, "estimateCovarianceFunctions", Config::OPTIONAL, "", ""))
    {
      estimateCovarianceFunctionVCE = TRUE;
      readConfig(config, "outputfileCovarianceFunctionPod", fileNameCovPod, Config::OPTIONAL, "", "covariance functions for along, cross, radial direction");
      endSequence(config);
    }
    if(readConfigSequence(config, "computeResiduals", Config::OPTIONAL, "", ""))
    {
      estimateResiduals = TRUE;
      readConfig(config, "outputfilePodResiduals", fileNameResidualsPod, Config::OPTIONAL, "", "");
      endSequence(config);
    }
    readConfig(config, "observation", observationMisc,        Config::MUSTSET,  "", "");
    if(readConfigSequence(config, "covariancePod", Config::MUSTSET, "", ""))
    {
      readConfig(config, "sigma",                        sigma0Pod,       Config::DEFAULT,  "1",  "apriori factor of covariance function");
      readConfig(config, "inputfileSigmasPerArc",        sigmaPodName,    Config::OPTIONAL, "",   "apriori different accuaries for each arc (multiplicated with sigma)");
      readConfig(config, "inputfileCovarianceFunction",  covPodName,      Config::OPTIONAL, "",   "approximate covariances in time");
      readConfig(config, "inputfileCovariancePodEpoch",  covPodEpochName, Config::OPTIONAL, "",   "3x3 epoch covariances");
      readConfig(config, "sampling",                     samplingPod,     Config::DEFAULT,  "30", "[seconds] sampling of the covariance function");
      endSequence(config);
    }
    readConfig(config, "inputfileArcList",    fileNameArcList,     Config::OPTIONAL, "",   "list to correspond points of time to arc numbers");
    readConfig(config, "adjustmentThreshold", adjustmentThreshold, Config::DEFAULT,  "0",  "Adjustment factor threshold: Iteration will be stopped once both SST and POD adjustment factors are under this threshold");
    readConfig(config, "iterationCount",      iterCount,           Config::DEFAULT,  "3",  "(maximum) number of iterations for the estimation of calibration parameter and error PSD");
    if(isCreateSchema(config)) return;

    // =============================================

    // init
    // ----
    const UInt arcCount        = observationMisc->arcCount();
    const UInt countAParameter = observationMisc->parameterCount();

    std::vector<UInt> arcsInterval = {0, arcCount};
    if(!fileNameArcList.empty())
    {
      logStatus<<"read arc list <"<<fileNameArcList<<">"<<Log::endl;
      std::vector<Time> timesInterval;
      readFileArcList(fileNameArcList, arcsInterval, timesInterval);
    }

    // init arc sigmas
    // ---------------
    sigmaPod = Vector(arcCount);
    for(UInt arcNo=0; arcNo<arcCount; arcNo++)
      sigmaPod(arcNo) = 1.0;
    if(!sigmaPodName.empty()) readFileMatrix(sigmaPodName, sigmaPod);
    if(sigmaPod.rows() != arcCount) throw(Exception("sigmasPerArc (POD) contains wrong number of arcs"));

    fileCovPodEpoch.open(covPodEpochName);

    // =============================================

    // Init normal equations
    // ---------------------
    logStatus<<"Init normal equations"<<Log::endl;
    std::vector<UInt> blockIndex(1, 0);
    blockIndex.push_back(countAParameter);
    normals.init(blockIndex);
    n = Matrix(normals.parameterCount(), observationMisc->rightSideCount());

    // ===================================================

    // parameter names
    // ---------------
    if(!fileNameParaName.empty() && Parallel::isMaster())
    {
      logStatus<<"write parameter names <"<<fileNameParaName<<">"<<Log::endl;
      observationMisc->setInterval(date2time(9999,1,1), date2time(9999,1,2));
      std::vector<ParameterName> paraNameStatic;
      observationMisc->parameterName(paraNameStatic);
      writeFileParameterName(fileNameParaName, paraNameStatic);
    }

    // =============================================

    // setup observation equations
    // ---------------------------
    logStatus<<"set up observation equations"<<Log::endl;
    observationArc.resize(arcCount);
    arcIndexN.resize(arcCount);
    arcIndexA.resize(arcCount);
    std::vector<UInt> processNo = Parallel::forEachInterval(arcCount, arcsInterval, [this](UInt arcNo) {computeObservationEquation(arcNo);});
    observationMisc = ObservationMiscPodPtr(nullptr);

    // =============================================

    // Determine max. length of ovariance functions
    // --------------------------------------------
    UInt covLengthPod = 0;
    for(UInt arcNo=0; arcNo<arcCount; arcNo++)
      if(observationArc.at(arcNo).times.size())
      {
        const UInt len = static_cast<UInt>(round((observationArc.at(arcNo).times.back()-observationArc.at(arcNo).times.at(0)).seconds()/samplingPod))+1;
        covLengthPod = std::max(covLengthPod, len);
      }
    Parallel::reduceMax(covLengthPod);
    Parallel::broadCast(covLengthPod);

    // =============================================

    // init covariance function
    // ------------------------
    CosTransformPod = Vce::cosTransform(covLengthPod);
    covFuncPod = Vce::readCovarianceFunction(covPodName, covLengthPod, 3, samplingPod);
    covFuncPod.column(1,3) *= pow(sigma0Pod,2);
    PsdPod = CosTransformPod * covFuncPod.column(1,3);

    // =============================================

    // Iteration
    // ---------
    for(UInt iter=0; iter<iterCount; iter++)
    {
      logInfo<<"starting iteration "<<iter+1<<Log::endl;
      bool thresholdReached = false;

      // solve normal equations
      // ----------------------
      if(countAParameter)
      {
        logStatus<<"accumulate system of normal equations"<<Log::endl;
        normals.setNull();
        n.setNull();
        lPl      = 0;
        obsCount = 0;
        Parallel::forEachProcess(arcCount, [this](UInt arcNo) {buildNormals(arcNo);}, processNo);

        // collect system of normal equations
        // ----------------------------------
        logStatus<<"collect system of normal equations"<<Log::endl;
        Parallel::reduceSum(n);
        Parallel::reduceSum(obsCount);
        Parallel::reduceSum(lPl);
        normals.reduceSum();

        // =============================================

        // Regularize not used parameters
        // ------------------------------
        for(UInt i=0; i<normals.blockCount(); i++)
          if(normals.isMyRank(i,i))
          {
            Matrix &N = normals.N(i,i);
            for(UInt k=0; k<N.rows(); k++)
              if(N(k,k) == 0)
                N(k,k) = 1.;
          }

        // =============================================

        // cholesky and forward step
        // -------------------------
        logStatus<<"solve system of normal equations"<<Log::endl;
        x = normals.solve(n, TRUE/*timing*/);
        Parallel::broadCast(x);
        if(Parallel::isMaster())
          logInfo<<"  aposteriori sigma = "<<sqrt((lPl-inner(x, n))/(obsCount-normals.parameterCount()))<<Log::endl;

        // N contains now the cholesky decomposition
        Wz = Vce::monteCarlo(normals.parameterCount(), 100);
        normals.triangularSolve(Wz);
        Parallel::broadCast(Wz);

        if(Parallel::isMaster() && !fileNameSolution.empty())
        {
          logStatus<<"write solution to <"<<fileNameSolution<<">"<<Log::endl;
          const UInt blockIndexStatic = 0;
          writeFileMatrix(fileNameSolution, x.row(normals.blockIndex(blockIndexStatic), normals.parameterCount()-normals.blockIndex(blockIndexStatic)));
        }

        if(!fileNameSigmax.empty())
        {
          logStatus<<"inverte cholesky matrix and write standard deviations to <"<<fileNameSigmax<<">"<<Log::endl;
          const UInt blockIndexStatic = 0;
          for(UInt i=blockIndexStatic; i<normals.blockCount(); i++)
            for(UInt k=i; k<normals.blockCount(); k++)
              if(normals.rank(i,k) != 0)
              {
                if(normals.isMyRank(i,k))
                  Parallel::send(normals.N(i,k), 0);
                else if(Parallel::isMaster())
                  Parallel::receive(normals.N(i,k), normals.rank(i,k));
              }
          if(Parallel::isMaster())
          {
            const UInt count = normals.parameterCount()-normals.blockIndex(blockIndexStatic);
            Matrix W(count, Matrix::TRIANGULAR);
            for(UInt i=blockIndexStatic; i<normals.blockCount(); i++)
              for(UInt k=i; k<normals.blockCount(); k++)
                copy(normals.N(i,k), W.slice(normals.blockIndex(i)-normals.blockIndex(blockIndexStatic), normals.blockIndex(k)-normals.blockIndex(blockIndexStatic), normals.blockSize(i), normals.blockSize(k)));
            inverse(W);
            Vector diagonal(count);
            for(UInt z=0; z<count; z++)
              diagonal(z) = sqrt(quadsum(W.slice(z,z,1,count-z)));
            writeFileMatrix(fileNameSigmax, diagonal);
          }
          for(UInt i=blockIndexStatic; i<normals.blockCount(); i++)
            for(UInt k=i; k<normals.blockCount(); k++)
              if(!normals.isMyRank(i,k))
                normals.N(i,k) = Matrix();
        } // if(!fileNameSigmax.empty())

      } // if(countAParameter)
      Parallel::barrier();

      // =============================================

      if(estimateResiduals)
      {
        logStatus<<"compute residuals"<<Log::endl;
        arcListPod.clear(); arcListPod.resize(arcCount, Epoch::ORBIT);
        Parallel::forEachProcess(arcCount,   [this](UInt arcNo) {computeResiduals(arcNo);},     processNo);
        Parallel::forEachProcess(arcListPod, [this](UInt arcNo) {return collectArcPod(arcNo);}, processNo);

        if(Parallel::isMaster() && (!fileNameResidualsPod.empty()))
        {
          logStatus<<"write residual file <"<<fileNameResidualsPod<<">"<<Log::endl;
          InstrumentFile::write(fileNameResidualsPod, arcListPod);
        }
      } // if(estimateResiduals)

      // =============================================

      if((estimateArcSigmas || estimateCovarianceFunctionVCE))
      {
        logStatus<<"compute redundancies"<<Log::endl;
        sigmaPodNew = Vector(arcCount);
        ePePod = redundancyPod = Matrix(covLengthPod, 3); // for x,y,z
        Parallel::forEachProcess(arcCount, [this](UInt arcNo) {computeRedundancies(arcNo);}, processNo);
      }

      // =============================================

      // sigmas per arc
      // --------------
      if(estimateArcSigmas)
      {
        Parallel::reduceSum(sigmaPodNew);
        if(Parallel::isMaster())
        {
          sigmaPod = sigmaPodNew;
          Double sigma0Pod = Vce::meanSigma(sigmaPod);
          sigmaPod *= 1./sigma0Pod;
          logInfo<<"  POD sigma per arc (median): "<<sigma0Pod<<Log::endl;

          if(!fileNameArcSigmaPod.empty())
          {
            logStatus<<"write arc sigma file <"<<fileNameArcSigmaPod<<">"<<Log::endl;
            writeFileMatrix(fileNameArcSigmaPod, sigmaPod);
          }
        }
        Parallel::broadCast(sigmaPod);
      } // if(estimateArcSigmas)

      // =============================================

      // estimate new PSD through variance component estimation
      // ------------------------------------------------------
      if(estimateCovarianceFunctionVCE)
      {
        Parallel::reduceSum(ePePod);
        Parallel::reduceSum(redundancyPod);

        logStatus<<"compute psd through variance component estimation"<<Log::endl;
        if(Parallel::isMaster())
        {
          Double maxFactor = 0;
          Vce::estimatePsd(ePePod, redundancyPod, PsdPod, maxFactor);
          logInfo<<"  max. PSD adjustment factor: "<<maxFactor<<Log::endl;
          if (maxFactor < adjustmentThreshold)
          {
            logStatus<<"  adjustment threshold "<<adjustmentThreshold<<" reached after iteration "<<iter+1<<"."<<Log::endl;
            thresholdReached = true;
          }
        } // if(Parallel::isMaster())
        Parallel::broadCast(PsdPod);
        Parallel::broadCast(thresholdReached);
        // compute new covariance function
        copy(CosTransformPod*PsdPod, covFuncPod.column(1,PsdPod.columns()));

        if(Parallel::isMaster() && !fileNameCovPod.empty())
        {
          logStatus<<"write covariance function file <"<<fileNameCovPod<<">"<<Log::endl;
          writeFileMatrix(fileNameCovPod, covFuncPod);
        }
      } // if(estimateCovarianceFunctionVCE)

      // =============================================

      // bail if the iteration threshold has been reached
      Parallel::broadCast(thresholdReached);
      if (thresholdReached)
        break;

      if((!estimateArcSigmas) && (!estimateCovarianceFunctionVCE))
        break; // iterations not needed
    } // for(iter)

  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void PreprocessingPod::computeObservationEquation(UInt arcNo)
{
  try
  {
    arcIndexN.at(arcNo).clear();
    arcIndexA.at(arcNo).clear();
    for(UInt idBlock=0; idBlock<normals.blockCount(); idBlock++)
    {
      arcIndexN.at(arcNo).push_back( idBlock );
      arcIndexA.at(arcNo).push_back( normals.blockIndex(idBlock) );
    }

    // compute observation equations
    // -----------------------------
    observationArc.at(arcNo) = observationMisc->computeArc(arcNo);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void PreprocessingPod::buildNormals(UInt arcNo)
{
  try
  {
    if(observationArc.at(arcNo).A.size() == 0)
      return;

    // Decorrelation
    // -------------
    Matrix Wl = observationArc.at(arcNo).l;
    Matrix WA = observationArc.at(arcNo).A;
    Matrix WB = observationArc.at(arcNo).B;
    Matrix W;
    ObservationSigmaArc sigmaEpoch;
    Covariance3dArc     covPodEpoch = fileCovPodEpoch.readArc(arcNo);
    CovariancePod::decorrelate(observationArc.at(arcNo).pod, sigmaPod(arcNo), sigmaEpoch, covPodEpoch, covFuncPod, W, {Wl, WA, WB});

    // eliminate arc dependent parameters
    // ----------------------------------
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
    // -------------
    this->lPl      += quadsum(l_bar);
    this->obsCount += l_bar.rows();
    Matrix n2 = A_bar.trans() * l_bar;
    Matrix N2(A_bar.columns(), Matrix::SYMMETRIC);
    rankKUpdate(1., A_bar, N2);
    fillSymmetric(N2);

    for(UInt i=0; i<arcIndexA.at(arcNo).size(); i++)
    {
      const UInt idxN1 = arcIndexN.at(arcNo).at(i);
      const UInt idxA1 = arcIndexA.at(arcNo).at(i);

      // right hand sides
      axpy(1., n2.row(idxA1, normals.blockSize(idxN1)), n.row(normals.blockIndex(idxN1), normals.blockSize(idxN1)));

      // normal matrix
      for(UInt k=0; k<arcIndexA.at(arcNo).size(); k++)
        if(arcIndexN.at(arcNo).at(k)>=idxN1)
        {
          const UInt idxN2 = arcIndexN.at(arcNo).at(k);
          const UInt idxA2 = arcIndexA.at(arcNo).at(k);
          if(normals.N(idxN1, idxN2).size() == 0)
            normals.N(idxN1, idxN2) = ((idxN1==idxN2) ? Matrix(normals.blockSize(idxN1), Matrix::SYMMETRIC)
                                                      : Matrix(normals.blockSize(idxN1), normals.blockSize(idxN2)));
          axpy(1., N2.slice(idxA1, idxA2, normals.blockSize(idxN1), normals.blockSize(idxN2)), normals.N(idxN1, idxN2));
        }
    } // for(i)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void PreprocessingPod::computeRedundancies(UInt arcNo)
{
  try
  {
    if(observationArc.at(arcNo).l.size() == 0)
      return;

    // count observations and calculate index
    // --------------------------------------
    const UInt countPod = observationArc.at(arcNo).times.size();

    // Decorrelation
    // -------------
    Matrix Wl = observationArc.at(arcNo).l;
    Matrix WA = observationArc.at(arcNo).A;
    Matrix WB = observationArc.at(arcNo).B;
    Matrix WPod;
    ObservationSigmaArc sigmaEpoch;
    Covariance3dArc     covPodEpoch = fileCovPodEpoch.readArc(arcNo);
    CovariancePod::decorrelate(observationArc.at(arcNo).pod, sigmaPod(arcNo), sigmaEpoch, covPodEpoch, covFuncPod, WPod, {Wl, WA, WB});

    // eliminate arc dependent parameters
    // ----------------------------------
    if(WB.size())
    {
      Vector tau = QR_decomposition(WB);
      QTransMult(WB, tau, Wl); // transform observations: l:= Q'l
      Wl.row(0, WB.columns()).setNull(); // residuals: remove WB*x
      QMult(WB, tau, Wl); // back transformation

      if(WA.size())
      {
        QTransMult(WB, tau, WA); // transform design matrix A:=Q'A
        WA.row(0, WB.columns()).setNull(); // residuals: remove WB*x
        QMult(WB, tau, WA); // back transformation
      }

      generateQ(WB, tau);
      WB.setType(Matrix::GENERAL);
    }

    // decorrelated residuals
    // ----------------------
    Matrix We = Wl;
    Matrix WAz(Wl.rows(), Wz.columns());
    for(UInt i=0; i<arcIndexA.at(arcNo).size(); i++)
    {
      const UInt idxN = arcIndexN.at(arcNo).at(i);
      const UInt idxA = arcIndexA.at(arcNo).at(i);
      matMult(-1., WA.column(idxA, normals.blockSize(idxN)),  x.row(normals.blockIndex(idxN), normals.blockSize(idxN)), We);
      matMult( 1., WA.column(idxA, normals.blockSize(idxN)), Wz.row(normals.blockIndex(idxN), normals.blockSize(idxN)), WAz);
    }

    // ============================================

    if(!estimateCovarianceFunctionVCE)
    {
      const Double redundancy = 3*countPod - quadsum(WAz) - quadsum(WB);
      sigmaPodNew(arcNo) = sqrt(quadsum(We)/redundancy) * sigmaPod(arcNo);

      return;
    } // if(!estimateCovarianceFunctions)

    // ============================================

    // Variance component estimation
    // -----------------------------
    std::vector<UInt> index(countPod);
    for(UInt i=0; i<index.size(); i++)
      index.at(i) = static_cast<UInt>(round((observationArc.at(arcNo).times.at(i)-observationArc.at(arcNo).times.at(0)).seconds()/samplingPod));
    Double ePeSum=0, redundancySum=0;
    Matrix R;
    Vector WWe;
    Vce::redundancy(WPod, We, WAz, WB, R, WWe);
    Vce::psd(R, WWe, index, sigmaPod(arcNo), CosTransformPod, PsdPod, ePePod, redundancyPod, ePeSum, redundancySum);
    sigmaPodNew(arcNo) = sqrt(ePeSum/redundancySum) * sigmaPod(arcNo);  // compute new sigma (for this arc)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void PreprocessingPod::computeResiduals(UInt arcNo)
{
  try
  {
    if(observationArc.at(arcNo).l.size() == 0)
      return;

    // Residuals
    // ---------
    Matrix e = observationArc.at(arcNo).l;
    for(UInt i=0; i<arcIndexA.at(arcNo).size(); i++)
    {
      const UInt idxN = arcIndexN.at(arcNo).at(i);
      const UInt idxA = arcIndexA.at(arcNo).at(i);
      matMult(-1., observationArc.at(arcNo).A.column(idxA, normals.blockSize(idxN)), x.row(normals.blockIndex(idxN), normals.blockSize(idxN)), e);
    }

    // eliminate arc dependent parameters
    // ----------------------------------
    if(observationArc.at(arcNo).B.size())
    {
      Matrix We = e;
      Matrix WB = observationArc.at(arcNo).B;

      Matrix WPod;
      ObservationSigmaArc sigmaEpoch;
      Covariance3dArc     covPodEpoch = fileCovPodEpoch.readArc(arcNo);
      CovariancePod::decorrelate(observationArc.at(arcNo).pod, sigmaPod(arcNo), sigmaEpoch, covPodEpoch, covFuncPod, WPod, {We, WB});

      Vector tau = QR_decomposition(WB);
      QTransMult(WB, tau, We); // transform observations: l:= Q'l
      Matrix y = We.row(0, tau.rows());
      triangularSolve(1., WB.row(0, tau.rows()), y);
      matMult(-1, observationArc.at(arcNo).B, y, e);
    }

    // create Pod arc
    // ---------------
    Arc arcPod;
    for(UInt i=0; i<observationArc.at(arcNo).times.size(); i++)
    {
      OrbitEpoch epoch;
      epoch.time     = observationArc.at(arcNo).times.at(i);
      epoch.position = Vector3d(e(3*i+0,0), e(3*i+1,0), e(3*i+2,0));
      arcListPod.at(arcNo).push_back(epoch);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
