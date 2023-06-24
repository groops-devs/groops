/***********************************************/
/**
* @file preprocessingDualSst.cpp
*
* @brief Estimate covariance function / arc weights for two simultaneous SST-ll observations.
*
* @author Torsten Mayer-Guerr
* @author Andreas Kvas
* @date 2019-11-25
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This programs processes satellite-to-satellite-tracking (SST) and orbit observations in a GRACE like configuration.
Four different observation groups are considered separately: two types of SST and POD1/POD2 for the two satellites.
This program works similar to \program{PreprocessingSst}, see there for details. Here only the settings explained,
which are different.

Both SST observation types are reduced by the same background models and the same impact
of accelerometer measurements. The covariance matrix of the reduced observations should not consider
the the instrument noise only (\configClass{covarianceSst1/2}{covarianceSstType}) but must
take the cross correlations \configClass{covarianceAcc}{covarianceSstType} into account.
The covariance matrix of the reduced observations is given by
\begin{equation}
  \M\Sigma(\begin{bmatrix} \Delta l_{SST1} \\ \Delta l_{SST2} \end{bmatrix})
  = \begin{bmatrix} \M\Sigma_{SST1} + \M\Sigma_{ACC} & \M\Sigma_{ACC} \\
                   \M\Sigma_{ACC} & \M\Sigma_{SST2} + \M\Sigma_{ACC}
    \end{bmatrix}.
\end{equation}
)";

/***********************************************/

#include "programs/program.h"
#include "base/polynomial.h"
#include "parallel/matrixDistributed.h"
#include "files/fileArcList.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"
#include "files/fileParameterName.h"
#include "classes/parametrizationGravity/parametrizationGravity.h"
#include "misc/observation/observationMiscDualSstVariational.h"
#include "misc/varianceComponentEstimation.h"
#include "misc/kalmanProcessing.h"
#include "misc/normalsShortTimeStaticLongTime.h"

/***** CLASS ***********************************/

/** @brief Estimate covariance function / arc weights.
* @ingroup programsGroup */
class PreprocessingDualSst
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);

private:
  ObservationMiscDualSstVariationalPtr observationMisc;
  std::vector<ObservationMiscDualSstVariational::Arc> observationArc;
  InstrumentFile                       fileCovEpochPod1, fileCovEpochPod2;
  std::vector<UInt>                    arcsInterval;
  std::vector<Time>                    timesInterval;

  static constexpr UInt TYPESIZE = 5;
  enum Type : UInt {SST1, SST2, ACC, POD1, POD2};
  std::array<std::string, TYPESIZE> typeName; // = {"SST1", "SST2, "ACC", "POD1", "POD2"};

  Bool estimateCovarianceFunctionVCE;
  Bool estimateArcSigmas;
  Bool estimateEpochSigmas;
  Bool estimateResiduals;
  Bool estimateSigmaShortTimeModel;

  // normal equations
  // ----------------
  NormalsShortTimeStaticLongTime normals;
  Matrix x;  // solution
  Matrix Wz; // monte carlo vector for redundancy computation

  // covariance
  std::array<UInt,   TYPESIZE> covColumns;
  std::array<Vector, TYPESIZE> sigma, sigmaNew;
  std::array<Matrix, TYPESIZE> covFunc, Psd, ePe, redundancy, CosTransform;
  std::array<Double, TYPESIZE> sampling;
  std::array<std::vector<ObservationSigmaArc>, TYPESIZE> arcListEpochSigma;
  std::array<std::vector<std::vector<Matrix>>, TYPESIZE> Cov; // Several independant matrices per arc
  std::array<Vector, TYPESIZE> sigmaCov;
  Double                       interpolationDegree;

  // residuals
  std::array<std::vector<Arc>, TYPESIZE> arcListResiduals;

  UInt findInterval(UInt arcNo) const;
  void decorrelate(UInt arcNo, std::array<Matrix, TYPESIZE> &Cov, std::array<Matrix, TYPESIZE> &W, const std::list<MatrixSlice> &A);
  void buildNormals(UInt arcNo);
  void computeRedundancies(UInt arcNo);
  void computeResiduals(UInt arcNo);
  void computeEpochSigmas(UInt arcNo);

  // Dual SST PSD VCE
  void varianceComponentEstimationPsd(const_MatrixSliceRef R, const_MatrixSliceRef WWe, const std::array<std::vector<Time>, TYPESIZE> &times,
                                      const std::array<std::vector<UInt>, TYPESIZE> &index, const std::array<Double, TYPESIZE> &sigma,
                                      const std::array<Matrix, TYPESIZE> &CosTransform, const std::array<Matrix, TYPESIZE> &Psd,
                                      std::array<Matrix, TYPESIZE> &ePe, std::array<Matrix, TYPESIZE> &redundancy,
                                      std::array<Double, TYPESIZE> &ePeSum, std::array<Double, TYPESIZE> &redundancySum);
};

GROOPS_REGISTER_PROGRAM(PreprocessingDualSst, PARALLEL, "Estimate covariance function / arc weights.", Preprocessing, Covariance, Residuals, KalmanFilter)

/***********************************************/
/***********************************************/

void PreprocessingDualSst::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    typeName   = {"SST1", "SST2", "ACC", "POD1", "POD2"};
    covColumns = {1, 1, 1, 3, 3};
    estimateCovarianceFunctionVCE = estimateArcSigmas = estimateEpochSigmas = estimateResiduals = FALSE;
    estimateSigmaShortTimeModel = FALSE;

    FileName                       fileNameSolution, fileNameSigmax, fileNameParaName;
    std::array<FileName, TYPESIZE> fileNameOutArcSigma, fileNameOutEpochSigma, fileNameOutCovFunc, fileNameResiduals;
    std::array<FileName, TYPESIZE> fileNameInArcSigma, fileNameInEpochSigma, fileNameInCovFunc;
    std::array<std::vector<FileName>, TYPESIZE> fileNamesInCov;
    std::array<FileName, TYPESIZE> fileNameInSigmasCov;
    FileName                       fileNameInCovEpochPod1,  fileNameInCovEpochPod2;
    std::array<Double, TYPESIZE>   sigma0; sigma0.fill(1.0);
    AutoregressiveModelSequencePtr arSequence;
    ParameterSelectorPtr           parameterShortTime;
    Double                         downweightPod;
    FileName                       fileNameArcList;
    UInt                           iterCount;
    std::string                    iterVariableName;
    UInt                           defaultBlockSize;

    readConfig(config, "outputfileSolution",      fileNameSolution, Config::OPTIONAL, "", "estimated parameter vector (static part only)");
    readConfig(config, "outputfileSigmax",        fileNameSigmax,   Config::OPTIONAL, "", "standard deviations of the parameters (sqrt of the diagonal of the inverse normal equation)");
    readConfig(config, "outputfileParameterName", fileNameParaName, Config::OPTIONAL, "", "estimated signal parameters (index is appended)");
    if(readConfigSequence(config, "estimateArcSigmas", Config::OPTIONAL, "", ""))
    {
      estimateArcSigmas = TRUE;
      readConfig(config, "outputfileSigmasPerArcSst1", fileNameOutArcSigma.at(SST1), Config::OPTIONAL, "", "accuracies of each arc (SST1)");
      readConfig(config, "outputfileSigmasPerArcSst2", fileNameOutArcSigma.at(SST2), Config::OPTIONAL, "", "accuracies of each arc (SST2)");
      readConfig(config, "outputfileSigmasPerArcAcc",  fileNameOutArcSigma.at(ACC),  Config::OPTIONAL, "", "accuracies of each arc (ACC)");
      readConfig(config, "outputfileSigmasPerArcPod1", fileNameOutArcSigma.at(POD1), Config::OPTIONAL, "", "accuracies of each arc (POD1)");
      readConfig(config, "outputfileSigmasPerArcPod2", fileNameOutArcSigma.at(POD2), Config::OPTIONAL, "", "accuracies of each arc (POD2)");
      endSequence(config);
    }
    if(readConfigSequence(config, "estimateEpochSigmas", Config::OPTIONAL, "", ""))
    {
      estimateEpochSigmas = TRUE;
      readConfig(config, "outputfileSigmasPerEpochSst1", fileNameOutEpochSigma.at(SST1), Config::OPTIONAL, "", "accuracies of each epoch (SST1)");
      readConfig(config, "outputfileSigmasPerEpochSst2", fileNameOutEpochSigma.at(SST2), Config::OPTIONAL, "", "accuracies of each epoch (SST2)");
      readConfig(config, "outputfileSigmasPerEpochAcc",  fileNameOutEpochSigma.at(ACC),  Config::OPTIONAL, "", "accuracies of each epoch (ACC)");
      readConfig(config, "outputfileSigmasPerEpochPod1", fileNameOutEpochSigma.at(POD1), Config::OPTIONAL, "", "accuracies of each epoch (POD1)");
      readConfig(config, "outputfileSigmasPerEpochPod2", fileNameOutEpochSigma.at(POD2), Config::OPTIONAL, "", "accuracies of each epoch (POD2)");
      endSequence(config);
    }
    if(readConfigSequence(config, "estimateCovarianceFunctions", Config::OPTIONAL, "", ""))
    {
      estimateCovarianceFunctionVCE = TRUE;
      readConfig(config, "outputfileCovarianceFunctionSst1", fileNameOutCovFunc.at(SST1), Config::OPTIONAL, "", "covariance function");
      readConfig(config, "outputfileCovarianceFunctionSst2", fileNameOutCovFunc.at(SST2), Config::OPTIONAL, "", "covariance function");
      readConfig(config, "outputfileCovarianceFunctionAcc",  fileNameOutCovFunc.at(ACC),  Config::OPTIONAL, "", "covariance function");
      readConfig(config, "outputfileCovarianceFunctionPod1", fileNameOutCovFunc.at(POD1), Config::OPTIONAL, "", "covariance functions for along, cross, radial direction");
      readConfig(config, "outputfileCovarianceFunctionPod2", fileNameOutCovFunc.at(POD2), Config::OPTIONAL, "", "covariance functions for along, cross, radial direction");
      endSequence(config);
    }
    if(readConfigSequence(config, "computeResiduals", Config::OPTIONAL, "", ""))
    {
      estimateResiduals = TRUE;
      readConfig(config, "outputfileSst1Residuals", fileNameResiduals.at(SST1), Config::OPTIONAL, "", "");
      readConfig(config, "outputfileSst2Residuals", fileNameResiduals.at(SST2), Config::OPTIONAL, "", "");
      readConfig(config, "outputfileAccResiduals",  fileNameResiduals.at(ACC),  Config::OPTIONAL, "", "");
      readConfig(config, "outputfilePod1Residuals", fileNameResiduals.at(POD1), Config::OPTIONAL, "", "");
      readConfig(config, "outputfilePod2Residuals", fileNameResiduals.at(POD2), Config::OPTIONAL, "", "");
      endSequence(config);
    }
    readConfig(config, "observation", observationMisc, Config::MUSTSET,  "", "");
    if(readConfigSequence(config, "covarianceSst1", Config::MUSTSET, "", ""))
    {
      readConfig(config, "sigma",                              sigma0.at(SST1),               Config::DEFAULT,  "1", "apriori factor of covariance function");
      readConfig(config, "inputfileSigmasPerArc",              fileNameInArcSigma.at(SST1),   Config::OPTIONAL, "",  "apriori different accuaries for each arc (multiplicated with sigma)");
      readConfig(config, "inputfileSigmasPerEpoch",            fileNameInEpochSigma.at(SST1), Config::OPTIONAL, "",  "apriori different accuaries for each epoch");
      readConfig(config, "inputfileCovarianceFunction",        fileNameInCovFunc.at(SST1),    Config::OPTIONAL, "",  "approximate covariances in time");
      readConfig(config, "inputfileCovarianceMatrixArc",       fileNamesInCov.at(SST1),       Config::OPTIONAL, "",  "Must be given per sst arc with correct dimensions.");
      readConfig(config, "inputfileSigmasCovarianceMatrixArc", fileNameInSigmasCov.at(SST1),  Config::OPTIONAL, "",  "Vector with one sigma for each <inputfileCovarianceMatrixArc>");
      readConfig(config, "sampling",                           sampling.at(SST1),             Config::DEFAULT,  "5", "[seconds] sampling of the covariance function");
      endSequence(config);
    }
    if(readConfigSequence(config, "covarianceSst2", Config::MUSTSET, "", ""))
    {
      readConfig(config, "sigma",                              sigma0.at(SST2),               Config::DEFAULT,  "1", "apriori factor of covariance function");
      readConfig(config, "inputfileSigmasPerArc",              fileNameInArcSigma.at(SST2),   Config::OPTIONAL, "",  "apriori different accuaries for each arc (multiplicated with sigma)");
      readConfig(config, "inputfileSigmasPerEpoch",            fileNameInEpochSigma.at(SST2), Config::OPTIONAL, "",  "apriori different accuaries for each epoch");
      readConfig(config, "inputfileCovarianceFunction",        fileNameInCovFunc.at(SST2),        Config::OPTIONAL, "",  "approximate covariances in time");
      readConfig(config, "inputfileCovarianceMatrixArc",       fileNamesInCov.at(SST2),       Config::OPTIONAL, "",  "Must be given per sst arc with correct dimensions.");
      readConfig(config, "inputfileSigmasCovarianceMatrixArc", fileNameInSigmasCov.at(SST2),  Config::OPTIONAL, "",  "Vector with one sigma for each <inputfileCovarianceMatrixArc>");
      readConfig(config, "sampling",                           sampling.at(SST2),             Config::DEFAULT,  "2", "[seconds] sampling of the covariance function");
      endSequence(config);
    }
    if(readConfigSequence(config, "covarianceAcc", Config::MUSTSET, "", ""))
    {
      readConfig(config, "sigma",                              sigma0.at(ACC),                Config::DEFAULT,  "1", "apriori factor of covariance function");
      readConfig(config, "inputfileSigmasPerArc",              fileNameInArcSigma.at(ACC),    Config::OPTIONAL, "",  "apriori different accuaries for each arc (multiplicated with sigma)");
      readConfig(config, "inputfileSigmasPerEpoch",            fileNameInEpochSigma.at(ACC),  Config::OPTIONAL, "",  "apriori different accuaries for each epoch");
      readConfig(config, "inputfileCovarianceFunction",        fileNameInCovFunc.at(ACC),     Config::OPTIONAL, "",  "approximate covariances in time");
      readConfig(config, "inputfileCovarianceMatrixArc",       fileNamesInCov.at(ACC),        Config::OPTIONAL, "",  "Must be given per sst arc with correct dimensions.");
      readConfig(config, "inputfileSigmasCovarianceMatrixArc", fileNameInSigmasCov.at(ACC),   Config::OPTIONAL, "",  "Vector with one sigma for each <inputfileCovarianceMatrixArc>");
      readConfig(config, "sampling",                           sampling.at(ACC),              Config::DEFAULT,  "5", "[seconds] sampling of the covariance function");
      endSequence(config);
    }
    if(readConfigSequence(config, "covariancePod1", Config::MUSTSET, "", ""))
    {
      readConfig(config, "sigma",                        sigma0.at(POD1),               Config::DEFAULT,  "1",  "apriori factor of covariance function");
      readConfig(config, "inputfileSigmasPerArc",        fileNameInArcSigma.at(POD1),   Config::OPTIONAL, "",   "apriori different accuaries for each arc (multiplicated with sigma)");
      readConfig(config, "inputfileSigmasPerEpoch",      fileNameInEpochSigma.at(POD1), Config::OPTIONAL, "",   "apriori different accuaries for each epoch");
      readConfig(config, "inputfileCovarianceFunction",  fileNameInCovFunc.at(POD1),    Config::OPTIONAL, "",   "approximate covariances in time");
      readConfig(config, "inputfileCovariancePodEpoch",  fileNameInCovEpochPod1,        Config::OPTIONAL, "",   "3x3 epoch covariances");
      readConfig(config, "sampling",                     sampling.at(POD1),             Config::DEFAULT,  "30", "[seconds] sampling of the covariance function");
      endSequence(config);
    }
    if(readConfigSequence(config, "covariancePod2", Config::MUSTSET, "", ""))
    {
      readConfig(config, "sigma",                        sigma0.at(POD2),               Config::DEFAULT,  "1",  "apriori factor of covariance function");
      readConfig(config, "inputfileSigmasPerArc",        fileNameInArcSigma.at(POD2),   Config::OPTIONAL, "",   "apriori different accuaries for each arc (multiplicated with sigma)");
      readConfig(config, "inputfileSigmasPerEpoch",      fileNameInEpochSigma.at(POD2), Config::OPTIONAL, "",   "apriori different accuaries for each epoch");
      readConfig(config, "inputfileCovarianceFunction",  fileNameInCovFunc.at(POD2),    Config::OPTIONAL, "",   "approximate covariances in time");
      readConfig(config, "inputfileCovariancePodEpoch",  fileNameInCovEpochPod2,        Config::OPTIONAL, "",   "3x3 epoch covariances");
      readConfig(config, "sampling",                     sampling.at(POD2),             Config::DEFAULT,  "30", "[seconds] sampling of the covariance function");
      endSequence(config);
    }
    if(readConfigSequence(config, "estimateShortTimeVariations", Config::OPTIONAL, "", "co-estimate short time gravity field variations"))
    {
      readConfig(config, "estimateSigma",               estimateSigmaShortTimeModel, Config::DEFAULT, "0", "estimate standard deviation via VCE");
      readConfig(config, "autoregressiveModelSequence", arSequence,                  Config::MUSTSET, "",  "AR model sequence for constraining short time gravity variations");
      readConfig(config, "parameterSelection",          parameterShortTime,          Config::MUSTSET, "",  "parameters describing the short time gravity field");
      endSequence(config);
    }
    readConfig(config, "downweightPod",          downweightPod,    Config::DEFAULT,  "1",    "downweight factor for POD");
    readConfig(config, "inputfileArcList",       fileNameArcList,  Config::OPTIONAL, "",     "list to correspond points of time to arc numbers");
    readConfig(config, "iterationCount",         iterCount,        Config::DEFAULT,  "3",    "(maximum) number of iterations for the estimation of calibration parameter and error PSD");
    readConfig(config, "variableNameIterations", iterVariableName, Config::OPTIONAL, "",     "All output fileNames in preprocessing iteration are expanded with this variable prior to writing to disk");
    readConfig(config, "defaultBlockSize",       defaultBlockSize, Config::DEFAULT,  "2048", "block size of static normal equation blocks");
    if(isCreateSchema(config)) return;

    // =============================================

    // init
    // ----
    interpolationDegree        = observationMisc->interpolationDegree;
    const UInt arcCount        = observationMisc->arcCount();
    const UInt countAParameter = observationMisc->parameterCount();

    arcsInterval  = {0, arcCount};
    timesInterval = {Time(), date2time(9999,1,1)};
    if(!fileNameArcList.empty())
    {
      logStatus<<"read arc list <"<<fileNameArcList<<">"<<Log::endl;
      readFileArcList(fileNameArcList, arcsInterval, timesInterval);
    }

    // init arc sigmas
    // ---------------
    for(UInt idType : {SST1, SST2, ACC, POD1, POD2})
    {
      sigma.at(idType) = Vector(arcCount, 1.0);
      if(!fileNameInArcSigma.at(idType).empty())
        readFileMatrix(fileNameInArcSigma.at(idType), sigma.at(idType));
      if(sigma.at(idType).rows()  != arcCount)
        throw(Exception("sigmasPerArc "+typeName.at(idType)+" contains wrong number of arcs"));
    }

    fileCovEpochPod1.open(fileNameInCovEpochPod1);
    fileCovEpochPod2.open(fileNameInCovEpochPod2);

    // ===================================================

    // normal equations of short time model
    // ------------------------------------
    UInt countShortTimeParameters = 0;
    std::vector<std::vector<std::vector<Matrix>>> normalsShortTime;
    if(arSequence)
    {
      logStatus<<"initialize normal equations of short time gravity field model"<<Log::endl;
      countShortTimeParameters = arSequence->dimension();
      normalsShortTime = arSequence->normalEquationSequence();
    }

    // init normal equations
    // ---------------------
    logStatus<<"initialize normal equations"<<Log::endl;
    normals.init(observationMisc, timesInterval, defaultBlockSize, comm, FALSE/*sortStateBeforeGravityParameter*/, countShortTimeParameters, parameterShortTime);

    // parameter names
    // ---------------
    if(!fileNameParaName.empty() && Parallel::isMaster(comm))
    {
      logStatus<<"write parameter names <"<<fileNameParaName<<">"<<Log::endl;
      writeFileParameterName(fileNameParaName, normals.parameterNames);
    }

    // =============================================

    // setup observation equations
    // ---------------------------
    logStatus<<"set up observation equations"<<Log::endl;
    observationArc.resize(arcCount);
    std::vector<UInt> processNo = Parallel::forEachInterval(arcCount, arcsInterval, [this](UInt arcNo)
    {
      const UInt idInterval = findInterval(arcNo);
      if(timesInterval.size())
        observationMisc->setInterval(timesInterval.at(idInterval), timesInterval.at(idInterval+1));
      observationArc.at(arcNo) = observationMisc->computeArc(arcNo);
    }, comm);
    observationMisc = nullptr;

    // =============================================

    // set used blocks
    // ---------------
    logStatus<<"setup normal equations"<<Log::endl;
    normals.setBlocks(arcsInterval);

    // =============================================

    // read SST covariance matrices
    // ----------------------------
    for(UInt idType : {SST1, SST2, ACC})
    {
      Cov.at(idType).resize(arcCount);
      if(fileNamesInCov.at(idType).size())
      {
        logStatus<<"read arc-wise covariance matrices"<<Log::endl;
        Parallel::forEachProcess(arcCount, [&](UInt arcNo)
        {
          Cov.at(idType).at(arcNo).resize(fileNamesInCov.at(idType).size());
          for(UInt i=0; i<fileNamesInCov.at(idType).size(); i++)
            readFileMatrix(fileNamesInCov.at(idType).at(i).appendBaseName(".arc"+arcNo%"%03i"s), Cov.at(idType).at(arcNo).at(i));
        }, processNo, comm);
      }

      sigmaCov.at(idType) = Vector(fileNamesInCov.at(idType).size(), 1.);
      if(!fileNameInSigmasCov.at(idType).empty())
        readFileMatrix(fileNameInSigmasCov.at(idType), sigmaCov.at(idType));
      if(sigmaCov.at(idType).rows() != fileNamesInCov.at(idType).size())
        throw(Exception("Number of sigmas not compatible with number of given arc-wise covariance matrices"));
    }

    // =============================================

    // Estimation of common component (ACC) only possible if both SSTs observed together
    // ---------------------------------------------------------------------------------
    UInt countArcBothSst = 0;
    for(UInt arcNo=0; arcNo<arcCount; arcNo++)
      if(observationArc.at(arcNo).times.at(SST1).size() && observationArc.at(arcNo).times.at(SST2).size())
        countArcBothSst++;
    Parallel::reduceSum(countArcBothSst, 0, comm);
    Parallel::broadCast(countArcBothSst, 0, comm);
    if(countArcBothSst == 0)
    {
      logWarningOnce<<"Estimation of common component (ACC) disabled: Need common SST1 and SST2 observations"<<Log::endl;
      for(UInt arcNo=0; arcNo<arcCount; arcNo++)
        observationArc.at(arcNo).times.at(ACC).clear();
    }

    // Initalize covariance functions
    // ------------------------------
    std::array<UInt, TYPESIZE> covLength; covLength.fill(0);
    for(UInt idType : {SST1, SST2, ACC, POD1, POD2})
    {
      // Determine max. length of ovariance functions
      for(UInt arcNo=0; arcNo<arcCount; arcNo++)
        if(observationArc.at(arcNo).times.at(idType).size())
          covLength.at(idType) = std::max(covLength.at(idType), static_cast<UInt>(std::round((observationArc.at(arcNo).times.at(idType).back()-observationArc.at(arcNo).times.at(idType).front()).seconds()/sampling.at(idType))+1));
      Parallel::reduceMax(covLength.at(idType), 0, comm);
      Parallel::broadCast(covLength.at(idType), 0, comm);

      // transformation PSD <-> covFunc
      CosTransform.at(idType) = Vce::cosTransform(covLength.at(idType));

      // init covariance function
      covFunc.at(idType) = Vce::readCovarianceFunction(fileNameInCovFunc.at(idType),  covLength.at(idType),  covColumns.at(idType), sampling.at(idType));
      covFunc.at(idType).column(1, covColumns.at(idType))  *= std::pow(sigma0.at(idType), 2);
      Psd.at(idType) = CosTransform.at(idType) * covFunc.at(idType).column(1, covColumns.at(idType));

      // init epoch sigmas
      arcListEpochSigma.at(idType).resize(arcCount);
      if(estimateEpochSigmas)
      {
        InstrumentFile file(fileNameInEpochSigma.at(idType));
        for(UInt arcNo=0; arcNo<arcCount; arcNo++)
          if(Parallel::myRank(comm) == processNo.at(arcNo))
          {
            arcListEpochSigma.at(idType).at(arcNo) = file.readArc(arcNo);
            if(!arcListEpochSigma.at(idType).at(arcNo).size())
              for(UInt i=0; i<observationArc.at(arcNo).times.at(idType).size(); i++)
              {
                ObservationSigmaEpoch epoch;
                epoch.time = observationArc.at(arcNo).times.at(idType).at(i);
                arcListEpochSigma.at(idType).at(arcNo).push_back(epoch);
              }
          }
      } // if(estimateEpochSigmas)
    } // for(idType)

    // =============================================

    // Iteration
    // ---------
    Double sigma2ShortTimeModel = 1.;
    for(UInt iter=0; iter<iterCount; iter++)
    {
      logInfo<<"starting iteration "<<iter<<Log::endl;
      VariableList variableIteration;
      if(iterVariableName.size())
        variableIteration.setVariable(iterVariableName, iter);

      // solve normal equations
      // ----------------------
      if(countAParameter)
      {
        logStatus<<"accumulate system of normal equations"<<Log::endl;
        normals.setNull();
        Parallel::forEachProcess(arcCount, [this](UInt arcNo) {buildNormals(arcNo);}, processNo, comm);
        logStatus<<"collect system of normal equations"<<Log::endl;
        normals.reduceSum();

        // add normals of short time model
        // -------------------------------
        if(normalsShortTime.size())
        {
          logStatus<<"add normals of short time model"<<Log::endl;
          normals.addShortTimeNormals(sigma2ShortTimeModel, normalsShortTime);
        }

        // cholesky and forward step
        // -------------------------
        logStatus<<"solve system of normal equations"<<Log::endl;
        const Double sigma = normals.solve(x, Wz);
        logInfo<<"  aposteriori sigma = "<<sigma<<Log::endl;

        if(Parallel::isMaster(comm) && !fileNameSolution.empty())
        {
          logStatus<<"write solution to <"<<fileNameSolution(variableIteration)<<">"<<Log::endl;
          writeFileMatrix(fileNameSolution(variableIteration), x);
        }

        if(!fileNameSigmax.empty())
        {
          logStatus<<"inverte cholesky matrix and write standard deviations to <"<<fileNameSigmax(variableIteration)<<">"<<Log::endl;
          Vector sigma = normals.parameterStandardDeviation();
          if(Parallel::isMaster(comm))
            writeFileMatrix(fileNameSigmax(variableIteration), sigma);
        } // if(!fileNameSigmax.empty())

        if(estimateSigmaShortTimeModel && normalsShortTime.size())
        {
          logStatus<<"compute standard deviation of short time gravity model"<<Log::endl;
          Double s2 = normals.estimateShortTimeNormalsVariance(sigma2ShortTimeModel, normalsShortTime, x, Wz);
          logInfo<<"  sigma: "<<std::sqrt(s2)<<Log::endl;
          if(!std::isnan(s2) && (s2 > 0))
            sigma2ShortTimeModel = s2;
        }
      } // if(countAParameter)
      Parallel::barrier(comm);

      // compute residuals
      // --------------------
      if(estimateResiduals || estimateEpochSigmas)
      {
        logStatus<<"compute residuals"<<Log::endl;
        for(UInt idType : {SST1, SST2, ACC, POD1, POD2})
        {
          arcListResiduals.at(idType).clear();
          arcListResiduals.at(idType).resize(arcCount);
        }
        Parallel::forEachProcess(arcCount, [this](UInt arcNo) {computeResiduals(arcNo);}, processNo, comm);
        for(UInt idType : {SST1, SST2, ACC, POD1, POD2})
        {
          Parallel::forEachProcess(arcListResiduals.at(idType), [this, idType](UInt arcNo) {return arcListResiduals.at(idType).at(arcNo);}, processNo, comm, FALSE/*timing*/);
          if(Parallel::isMaster(comm) && (!fileNameResiduals.at(idType).empty()))
          {
            logStatus<<"write residual file <"<<fileNameResiduals.at(idType)(variableIteration)<<">"<<Log::endl;
            InstrumentFile::write(fileNameResiduals.at(idType)(variableIteration), arcListResiduals.at(idType));
          }
        }
      }

      // compute redundancies
      // --------------------
      if(estimateArcSigmas || estimateCovarianceFunctionVCE)
      {
        logStatus<<"compute redundancies"<<Log::endl;
        for(UInt idType : {SST1, SST2, ACC, POD1, POD2})
        {
          sigmaNew.at(idType) = Vector(arcCount);
          ePe.at(idType) = redundancy.at(idType) = Matrix(covLength.at(idType), covColumns.at(idType));
        }
        Parallel::forEachProcess(arcCount, [this](UInt arcNo) {computeRedundancies(arcNo);}, processNo, comm);
      }

      // sigmas per arc
      // --------------
      if(estimateArcSigmas)
      {
        for(UInt idType : {SST1, SST2, ACC, POD1, POD2})
        {
          Parallel::reduceSum(sigmaNew.at(idType), 0, comm);
          if(Parallel::isMaster(comm))
          {
            sigma.at(idType) = (1./Vce::meanSigma(sigmaNew.at(idType))) * sigmaNew.at(idType);
            logInfo<<"  "<<typeName.at(idType)<<" sigma per arc (median): "<<Vce::meanSigma(sigmaNew.at(idType))<<Log::endl;
            for(UInt arcNo=0; arcNo<sigma.at(idType).size(); arcNo++)
              if(sigma.at(idType)(arcNo) < 0)
                sigma.at(idType)(arcNo) = 1;
          }
          Parallel::broadCast(sigma.at(idType), 0, comm);
        }

        for(UInt idType : {SST1, SST2, ACC, POD1, POD2})
          if(Parallel::isMaster(comm) && !fileNameOutArcSigma.at(idType).empty())
          {
            logStatus<<"write arc sigma file <"<<fileNameOutArcSigma.at(idType)(variableIteration)<<">"<<Log::endl;
            writeFileMatrix(fileNameOutArcSigma.at(idType)(variableIteration), sigma.at(idType));
          }
      }

      // sigmas per epoch
      // --------------
      if(estimateEpochSigmas)
      {
        logStatus<<"estimate epoch sigmas"<<Log::endl;
        Parallel::forEachProcess(arcCount, [this](UInt arcNo) {computeEpochSigmas(arcNo);}, processNo, comm);
        for(UInt idType : {SST1, SST2, ACC, POD1, POD2})
        {
          Parallel::forEachProcess(arcListEpochSigma.at(idType), [this, idType](UInt arcNo) {return arcListEpochSigma.at(idType).at(arcNo);}, processNo, comm, FALSE/*timing*/);
          if(Parallel::isMaster(comm) && !fileNameOutEpochSigma.at(idType).empty())
          {
            logStatus<<"write epoch sigma file <"<<fileNameOutEpochSigma.at(idType)(variableIteration)<<">"<<Log::endl;
            InstrumentFile::write(fileNameOutEpochSigma.at(idType)(variableIteration), arcListEpochSigma.at(idType));
          }
        }
        Parallel::barrier(comm);
      } // if(estimateEpochSigmas)

      // estimate new PSD through variance component estimation
      // ------------------------------------------------------
      if(estimateCovarianceFunctionVCE)
      {
        logStatus<<"estimate PSDs"<<Log::endl;
        for(UInt idType : {SST1, SST2, ACC, POD1, POD2})
        {
          Parallel::reduceSum(ePe.at(idType), 0, comm);
          Parallel::reduceSum(redundancy.at(idType), 0, comm);
          if(Parallel::isMaster(comm))
          {
            Double maxFactor = 0;
            Vce::estimatePsd(ePe.at(idType), redundancy.at(idType), Psd.at(idType), maxFactor);
            if((idType == POD1) || (idType == POD2))
              maxFactor /= downweightPod;
            logInfo<<"  max. PSD adjustment factor "<<typeName.at(idType)<<": "<<maxFactor<<Log::endl;
          }
          Parallel::broadCast(Psd.at(idType), 0, comm);
          // compute new covariance function
          copy(CosTransform.at(idType) * Psd.at(idType), covFunc.at(idType).column(1, covColumns.at(idType)));
        }

        // Write covariance function to file
        for(UInt idType : {SST1, SST2, ACC, POD1, POD2})
          if(Parallel::isMaster(comm) && !fileNameOutCovFunc.at(idType).empty())
          {
            logStatus<<"write covariance function file <"<<fileNameOutCovFunc.at(idType)(variableIteration)<<">"<<Log::endl;
            writeFileMatrix(fileNameOutCovFunc.at(idType)(variableIteration), covFunc.at(idType));
          }
      }

      // downweight POD
      // --------------
      for(UInt idType : {POD1, POD2})
      {
        Psd.at(idType) *= std::pow(downweightPod, 2);
        covFunc.at(idType).column(1, covColumns.at(idType)) *= std::pow(downweightPod, 2);
      }
    } // for(iter)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

UInt PreprocessingDualSst::findInterval(UInt arcNo) const
{
  for(UInt idInterval=0; idInterval+1<arcsInterval.size(); idInterval++)
    if(arcNo < arcsInterval.at(idInterval+1))
      return idInterval;
  return 0;
}

/***********************************************/

void PreprocessingDualSst::decorrelate(UInt arcNo, std::array<Matrix, TYPESIZE> &Cov, std::array<Matrix, TYPESIZE> &W, const std::list<MatrixSlice> &A)
{
  try
  {
    const auto &times = observationArc.at(arcNo).times;

    // count observations and calculate index
    // --------------------------------------
    std::array<UInt, TYPESIZE> idx;
    std::array<std::list<MatrixSlice>, TYPESIZE> WA; // type specific slices of A
    UInt obsCount = 0;
    for(UInt idType : {SST1, SST2, POD1, POD2})
    {
      idx.at(idType) = obsCount;
      obsCount      += covColumns.at(idType) * times.at(idType).size();
      for(const MatrixSlice &a : A)
        WA.at(idType).push_back(a.row(idx.at(idType), covColumns.at(idType) * times.at(idType).size()));
    }

    if(times.at(SST1).size() || times.at(SST2).size())
    {
      std::list<MatrixSlice> WA;
      for(const MatrixSlice &a : A)
        WA.push_back(a.row(idx.at(SST1), times.at(SST1).size()+times.at(SST2).size()));

      for(UInt idType : {SST1, SST2, ACC})
        if(this->Cov.at(idType).at(arcNo).size())
        {
          Cov.at(idType) = std::pow(sigmaCov.at(idType).at(0), 2) * this->Cov.at(idType).at(arcNo).at(0);
          for(UInt i=1; i<this->Cov.at(idType).at(arcNo).size(); i++)
            axpy(std::pow(sigmaCov.at(idType).at(i), 2), this->Cov.at(idType).at(arcNo).at(i), Cov.at(idType));
        }

      for(UInt idType : {SST1, SST2, ACC})
        CovarianceSst::covariance(times.at(idType), sigma.at(idType)(arcNo), arcListEpochSigma.at(idType).at(arcNo), covFunc.at(idType), Cov.at(idType));
      W.at(ACC) = ObservationMiscDualSstVariational::decorrelate(times.at(SST1), times.at(SST2), times.at(ACC), Cov.at(SST1), Cov.at(SST2), Cov.at(ACC), interpolationDegree, WA);
    }

    if(times.at(POD1).size())
      W.at(POD1) = CovariancePod::decorrelate(observationArc.at(arcNo).pod1, sigma.at(POD1)(arcNo), arcListEpochSigma.at(POD1).at(arcNo), fileCovEpochPod1.readArc(arcNo), covFunc.at(POD1), WA.at(POD1));

    if(times.at(POD2).size())
      W.at(POD2) = CovariancePod::decorrelate(observationArc.at(arcNo).pod2, sigma.at(POD2)(arcNo), arcListEpochSigma.at(POD2).at(arcNo), fileCovEpochPod2.readArc(arcNo), covFunc.at(POD2), WA.at(POD2));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void PreprocessingDualSst::buildNormals(UInt arcNo)
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
    std::array<Matrix, TYPESIZE> Cov, W;
    decorrelate(arcNo, Cov, W, {Wl, WA, WB});

    normals.accumulate(findInterval(arcNo), Wl, WA, WB);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void PreprocessingDualSst::computeRedundancies(UInt arcNo)
{
  try
  {
    if(observationArc.at(arcNo).l.size() == 0)
      return;

    // count observations and calculate index
    // --------------------------------------
    const auto &times = observationArc.at(arcNo).times;
    std::array<UInt, TYPESIZE> count, idx;
    UInt obsCount = 0;
    for(UInt idType : {SST1, SST2, POD1, POD2}) // without ACC
    {
      count.at(idType) = times.at(idType).size();
      idx.at(idType)   = obsCount;
      obsCount        += covColumns.at(idType) * count.at(idType);
    }
    count.at(ACC) = times.at(ACC).size();

    // Decorrelation
    // -------------
    Matrix Wl = observationArc.at(arcNo).l;
    Matrix WA = observationArc.at(arcNo).A;
    Matrix WB = observationArc.at(arcNo).B;
    std::array<Matrix, TYPESIZE> Cov, W;
    decorrelate(arcNo, Cov, W, {Wl, WA, WB});

    // eliminate arc dependent parameters
    // ----------------------------------
    if(WB.size())
    {
      Vector tau = QR_decomposition(WB);
      QTransMult(WB, tau, Wl);             // transform observations: l:= Q'l
      Wl.row(0, WB.columns()).setNull();   // residuals: remove WB*x
      QMult(WB, tau, Wl);                  // back transformation
      if(WA.size())
      {
        QTransMult(WB, tau, WA);           // transform design matrix A:=Q'A
        WA.row(0, WB.columns()).setNull(); // residuals: remove WB*x
        QMult(WB, tau, WA);                // back transformation
      }
      generateQ(WB, tau);
    }

    // decorrelated residuals
    // ----------------------
    Matrix We = Wl;
    Matrix WAz(Wl.rows(), Wz.columns());
    const UInt idInterval = findInterval(arcNo);
    normals.designMatMult(idInterval, -1., WA, x,  We);
    normals.designMatMult(idInterval, +1., WA, Wz, WAz);

    // ============================================

    // variance component estimation
    // -----------------------------
    Matrix R;
    Vector WWe;
    std::array<std::vector<UInt>, TYPESIZE> index;
    std::array<Double, TYPESIZE> ePeSum, redundancySum; ePeSum.fill(0); redundancySum.fill(0);
    for(UInt idType : {SST1, SST2, ACC, POD1, POD2})
    {
      index.at(idType).resize(times.at(idType).size());
      for(UInt i=0; i<index.at(idType).size(); i++)
        index.at(idType).at(i) = static_cast<UInt>(std::round((times.at(idType).at(i)-times.at(idType).front()).seconds()/sampling.at(idType)));
    }

    if(count.at(SST1) || count.at(SST2))
    {
      Vce::redundancy(W.at(ACC), We.row(idx.at(SST1), count.at(SST1)+count.at(SST2)),
                      WAz.row(idx.at(SST1), count.at(SST1)+count.at(SST2)), WB.row(idx.at(SST1), count.at(SST1)+count.at(SST2)), R, WWe);
      varianceComponentEstimationPsd(R, WWe, times, index, {sigma.at(SST1)(arcNo), sigma.at(SST2)(arcNo), sigma.at(ACC)(arcNo), 0., 0.},
                                     CosTransform, Psd, ePe, redundancy, ePeSum, redundancySum);
    }

    for(UInt idType : {POD1, POD2})
      if(count.at(idType))
      {
        Vce::redundancy(W.at(idType), We.row(idx.at(idType), 3*count.at(idType)),
                        WAz.row(idx.at(idType), 3*count.at(idType)), WB.row(idx.at(idType), 3*count.at(idType)), R, WWe);
        Vce::psd(R, WWe, index.at(idType), sigma.at(idType)(arcNo), CosTransform.at(idType), Psd.at(idType),
                 ePe.at(idType), redundancy.at(idType), ePeSum.at(idType), redundancySum.at(idType));
      }

    for(UInt idType : {SST1, SST2, ACC, POD1, POD2}) // compute new sigma (for this arc)
      if(redundancySum.at(idType) > 0.3)
        sigmaNew.at(idType)(arcNo) = std::sqrt(ePeSum.at(idType)/redundancySum.at(idType)) * sigma.at(idType)(arcNo);

    // estimation of ACC arc sigmas not possible if not both SSTs are observed
    if((!count.at(SST1) || !count.at(SST2)) && sigmaNew.at(ACC)(arcNo))
      sigmaNew.at(ACC)(arcNo) = -1;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void PreprocessingDualSst::varianceComponentEstimationPsd(const_MatrixSliceRef R, const_MatrixSliceRef WWe,
                                                          const std::array<std::vector<Time>, TYPESIZE> &times,
                                                          const std::array<std::vector<UInt>, TYPESIZE> &index, const std::array<Double, TYPESIZE> &sigma,
                                                          const std::array<Matrix, TYPESIZE> &CosTransform, const std::array<Matrix, TYPESIZE> &Psd,
                                                          std::array<Matrix, TYPESIZE> &ePe, std::array<Matrix, TYPESIZE> &redundancy,
                                                          std::array<Double, TYPESIZE> &ePeSum, std::array<Double, TYPESIZE> &redundancySum)
{
  try
  {
    std::array<Vector, TYPESIZE> e, r;
    std::array<UInt, 2> idx = {0, index.at(SST1).size()};
    for(UInt idType : {SST1, SST2})
      if(CosTransform.at(idType).rows())
      {
        e.at(idType) = r.at(idType) = Vector(CosTransform.at(idType).rows());
        for(UInt i=0; i<index.at(idType).size(); i++)
          for(UInt k=i; k<index.at(idType).size(); k++)
          {
            e.at(idType)(index.at(idType).at(k)-index.at(idType).at(i)) += WWe(i+idx.at(idType),0) * WWe(k+idx.at(idType),0);
            r.at(idType)(index.at(idType).at(k)-index.at(idType).at(i)) += R(i+idx.at(idType), k+idx.at(idType));
          }
        e.at(idType).row(1, e.at(idType).rows()-1) *= 2.; // consider lower triangular of matrix
        r.at(idType).row(1, r.at(idType).rows()-1) *= 2.;
      }

    // mixed term: ACC
    // ---------------
    if(CosTransform.at(ACC).rows())
    {
      e.at(ACC) = r.at(ACC) = Vector(CosTransform.at(ACC).rows());
      Polynomial polynomial(times.at(ACC), interpolationDegree);
      fillSymmetric(R);
      Matrix WWe1 = polynomial.interpolate(times.at(SST1), WWe.row(idx.at(SST1), index.at(SST1).size()),        1, 0, TRUE/*adjoint*/);  // F1'*WWe1
      Matrix WWe2 = polynomial.interpolate(times.at(SST2), WWe.row(idx.at(SST2), index.at(SST2).size()),        1, 0, TRUE/*adjoint*/);  // F2'*WWe2
      Matrix R1   = polynomial.interpolate(times.at(SST1), R.row(idx.at(SST1), index.at(SST1).size()),          1, 0, TRUE/*adjoint*/);  // F1'*(R11,R12)
      Matrix R12  = polynomial.interpolate(times.at(SST2), R1.trans().row(idx.at(SST2), index.at(SST2).size()), 1, 0, TRUE/*adjoint*/);  // F2'*R12'*F1
             R1   = polynomial.interpolate(times.at(SST1), R1.trans().row(idx.at(SST1), index.at(SST1).size()), 1, 0, TRUE/*adjoint*/);  // F1'*R11*F1
      Matrix R2   = polynomial.interpolate(times.at(SST2), R.slice(idx.at(SST2), idx.at(SST2), index.at(SST2).size(), index.at(SST2).size()), 1, 0, TRUE); // F2'*R22
             R2   = polynomial.interpolate(times.at(SST2), R2.trans(),                                          1, 0, TRUE/*adjoint*/);  // F2'*R22*F2

      for(UInt i=0; i<index.at(ACC).size(); i++)
        for(UInt k=0; k<index.at(ACC).size(); k++)
        {
          const UInt idiff = (i < k) ? (index.at(ACC).at(k)-index.at(ACC).at(i)) : (index.at(ACC).at(i)-index.at(ACC).at(k));
          e.at(ACC)(idiff) += WWe1(i,0)*WWe1(k,0) +  WWe2(i,0)*WWe2(k,0) + WWe1(i,0)*WWe2(k,0) + WWe2(i,0)*WWe1(k,0);
          r.at(ACC)(idiff) += R1(i,k) + R2(i,k) + R12(i,k) + R12(k,i);
        }
    }

    for(UInt idType : {SST1, SST2, ACC})
      for(UInt idFreq=0; idFreq<Psd.at(idType).rows(); idFreq++)
      {
        const Double ePeTmp        = std::pow(sigma.at(idType), 2) * Psd.at(idType)(idFreq, 0) * inner(e.at(idType), CosTransform.at(idType).column(idFreq));
        const Double redundancyTmp = std::pow(sigma.at(idType), 2) * Psd.at(idType)(idFreq, 0) * inner(r.at(idType), CosTransform.at(idType).column(idFreq));
        ePe.at(idType)(idFreq, 0)        += ePeTmp;
        redundancy.at(idType)(idFreq, 0) += redundancyTmp;
        ePeSum.at(idType)                += ePeTmp;
        redundancySum.at(idType)         += redundancyTmp;
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void PreprocessingDualSst::computeResiduals(UInt arcNo)
{
  try
  {
    if(observationArc.at(arcNo).l.size() == 0)
      return;

    // count observations and calculate index
    // --------------------------------------
    const auto &times = observationArc.at(arcNo).times;
    std::array<UInt, TYPESIZE> count, idx;
    UInt obsCount = 0;
    for(UInt idType : {SST1, SST2, POD1, POD2}) // without ACC
    {
      count.at(idType) = times.at(idType).size();
      idx.at(idType)   = obsCount;
      obsCount        += covColumns.at(idType) * count.at(idType);
    }
    count.at(ACC) = times.at(ACC).size();

    // Residuals
    // ---------
    Matrix e = observationArc.at(arcNo).l;
    normals.designMatMult(findInterval(arcNo), -1., observationArc.at(arcNo).A, x, e);

    // eliminate arc dependent parameters
    // ----------------------------------
    std::array<Matrix, TYPESIZE> Cov, W;
    if(observationArc.at(arcNo).B.size())
    {
      Matrix We = e;
      Matrix WB = observationArc.at(arcNo).B;
      decorrelate(arcNo, Cov, W, {We, WB});
      Vector tau = QR_decomposition(WB);
      QTransMult(WB, tau, We); // transform observations: l:= Q'l
      Matrix y = We.row(0, tau.rows());
      triangularSolve(1., WB.row(0, tau.rows()), y);
      matMult(-1, observationArc.at(arcNo).B, y, e);
    }
    else
      decorrelate(arcNo, Cov, W, {});

    // residuals
    // ---------
    std::array<Matrix, TYPESIZE> residuals;
    if(Cov.at(ACC).size())
    {
      Polynomial polynomial(times.at(ACC), interpolationDegree);
      Matrix WWe = e.row(idx.at(SST1), count.at(SST1)+count.at(SST2));
      triangularSolve(1., W.at(ACC).trans(), WWe);
      triangularSolve(1., W.at(ACC),         WWe);
      residuals.at(SST1) = Cov.at(SST1) * WWe.row(idx.at(SST1), count.at(SST1)); // collocation for residual separation
      residuals.at(SST2) = Cov.at(SST2) * WWe.row(idx.at(SST2), count.at(SST2));
      residuals.at(ACC)  = Cov.at(ACC)  * (polynomial.interpolate(times.at(SST1), WWe.row(idx.at(SST1), count.at(SST1)), 1, 0, TRUE/*adjoint*/) +
                                          polynomial.interpolate(times.at(SST2), WWe.row(idx.at(SST2), count.at(SST2)), 1, 0, TRUE/*adjoint*/));
    }
    else
    {
      residuals.at(SST1) = e.row(idx.at(SST1), count.at(SST1)); // collocation for residual separation
      residuals.at(SST2) = e.row(idx.at(SST2), count.at(SST2));
    }
    residuals.at(POD1) = e.row(idx.at(POD1), covColumns.at(POD1) * count.at(POD1));
    residuals.at(POD2) = e.row(idx.at(POD2), covColumns.at(POD2) * count.at(POD2));

    // create Sst arc
    // --------------
    for(UInt idType : {SST1, SST2, ACC})
      for(UInt i=0; i<count.at(idType); i++)
      {
        SatelliteTrackingEpoch epoch;
        epoch.time  = times.at(idType).at(i);
        epoch.range = epoch.rangeRate = epoch.rangeAcceleration = residuals.at(idType)(i,0);
        arcListResiduals.at(idType).at(arcNo).push_back(epoch);
      }

    // create Pod arcs
    // ---------------
    for(UInt idType : {POD1, POD2})
      for(UInt i=0; i<count.at(idType); i++)
      {
        OrbitEpoch epoch;
        epoch.time     = times.at(idType).at(i);
        epoch.position = Vector3d(residuals.at(idType).row(3*i, 3));
        arcListResiduals.at(idType).at(arcNo).push_back(epoch);
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void PreprocessingDualSst::computeEpochSigmas(UInt arcNo)
{
  try
  {
    const Double huber = 2.5;

    for(UInt idType : {SST1, SST2, ACC})
      if(covFunc.at(idType).size())
      {
        const Double threshold2 = std::pow(huber*sigma.at(idType)(arcNo), 2) * covFunc.at(idType)(0,1);
        for(UInt i=0; i<observationArc.at(arcNo).times.at(idType).size(); i++)
        {
          const Double e2 = std::pow(dynamic_cast<const SatelliteTrackingEpoch &>(arcListResiduals.at(idType).at(arcNo).at(i)).rangeRate, 2);
          arcListEpochSigma.at(idType).at(arcNo).at(i).sigma = 0.;
          if(e2 > threshold2)
            arcListEpochSigma.at(idType).at(arcNo).at(i).sigma = std::sqrt(e2-threshold2);
        }
      }

    for(UInt idType : {POD1, POD2})
      if(covFunc.at(idType).size())
      {
        const Double threshold2 = std::pow(huber*sigma.at(idType)(arcNo), 2) * sum(covFunc.at(idType).slice(0, 1, 1, 3));
        for(UInt i=0; i<observationArc.at(arcNo).times.at(idType).size(); i++)
        {
          const Double e2 = dynamic_cast<const OrbitEpoch &>(arcListResiduals.at(idType).at(arcNo).at(i)).position.quadsum();
          arcListEpochSigma.at(idType).at(arcNo).at(i).sigma = 0.;
          if(e2 > threshold2)
            arcListEpochSigma.at(idType).at(arcNo).at(i).sigma = std::sqrt((e2-threshold2)/3);
        }
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
