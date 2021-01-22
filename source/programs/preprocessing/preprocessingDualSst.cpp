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
Four different observation groups are considered separatly: two types of SST and POD1/POD2 for the two satellites.
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
#include "base/fourier.h"
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
  InstrumentFile    fileCovEpochPod1, fileCovEpochPod2;
  std::vector<UInt> arcsInterval;
  std::vector<Time> timesInterval;

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
  // ----------
  Vector  sigmaSst1,    sigmaSst2,    sigmaAcc,    sigmaPod1,    sigmaPod2;
  Vector  sigmaSst1New, sigmaSst2New, sigmaAccNew, sigmaPod1New, sigmaPod2New;
  std::vector<ObservationSigmaArc> arcListEpochSigmaSst1, arcListEpochSigmaSst2, arcListEpochSigmaAcc, arcListEpochSigmaPod1, arcListEpochSigmaPod2;
  std::vector<std::vector<Matrix>> CovSst1, CovSst2, CovAcc;
  Vector  sigmaCovSst1, sigmaCovSst2, sigmaCovAcc;
  Vector  covMatrixSst1Sigmas, ePeCovMatrixSst1;
  Vector  covMatrixSst2Sigmas, ePeCovMatrixSst2;
  Vector  covMatrixAccSigmas,  ePeCovMatrixAcc;
  Matrix  covFuncSst1, covFuncSst2, covFuncAcc, covFuncPod1, covFuncPod2;
  Matrix  PsdSst1, PsdSst2, PsdAcc, PsdPod1, PsdPod2;
  Matrix  ePeSst1, ePeSst2, ePeAcc, ePePod1, ePePod2;
  Matrix  redundancySst1, redundancySst2, redundancyAcc, redundancyPod1, redundancyPod2; // one row for each frequency, one column for each component
  Matrix  CosTransformSst, CosTransformPod1, CosTransformPod2;
  Double  samplingSst, samplingPod1, samplingPod2;

  // residuals
  std::vector<SatelliteTrackingArc> arcListResidualsSst1, arcListResidualsSst2, arcListResidualsAcc;
  std::vector<OrbitArc> arcListResidualsPod1, arcListResidualsPod2;

  UInt findInterval(UInt arcNo) const;
  void decorrelate(UInt arcNo, UInt countSst, UInt countPod1, UInt countPod2, Matrix &CovSst1, Matrix &CovSst2, Matrix &CovAcc,
                   Matrix &WSst, Matrix &WPod1, Matrix &WPod2, const std::list<MatrixSlice> &A);
  void buildNormals(UInt arcNo);
  void computeRedundancies(UInt arcNo);
  void computeResiduals(UInt arcNo);
  void computeEpochSigmas(UInt arcNo);

  // Dual SST PSD VCE
  void varianceComponentEstimationPsd(const_MatrixSliceRef R, const_MatrixSliceRef WWe, const_MatrixSliceRef CosTransform,
                                      Double sigma1, Double sigma2, Double sigma12,
                                      const_MatrixSliceRef Psd1, const_MatrixSliceRef Psd2, const_MatrixSliceRef Psd12,
                                      MatrixSliceRef ePe1,  MatrixSliceRef redundancy1,  Double &ePe1Sum,  Double &redundancy1Sum,
                                      MatrixSliceRef ePe2,  MatrixSliceRef redundancy2,  Double &ePe2Sum,  Double &redundancy2Sum,
                                      MatrixSliceRef ePe12, MatrixSliceRef redundancy12, Double &ePe12Sum, Double &redundancy12Sum);
};

GROOPS_REGISTER_PROGRAM(PreprocessingDualSst, PARALLEL, "Estimate covariance function / arc weights.", Preprocessing, Covariance, Residuals, KalmanFilter)

/***********************************************/
/***********************************************/

void PreprocessingDualSst::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    estimateCovarianceFunctionVCE = estimateArcSigmas = estimateEpochSigmas = estimateResiduals = FALSE;
    estimateSigmaShortTimeModel = FALSE;

    FileName fileNameSolution,          fileNameSigmax,            fileNameParaName;
    FileName fileNameOutArcSigmaSst1,   fileNameOutArcSigmaSst2,   fileNameOutArcSigmaAcc,   fileNameOutArcSigmaPod1,   fileNameOutArcSigmaPod2;
    FileName fileNameOutEpochSigmaSst1, fileNameOutEpochSigmaSst2, fileNameOutEpochSigmaAcc, fileNameOutEpochSigmaPod1, fileNameOutEpochSigmaPod2;
    FileName fileNameOutCovSst1,        fileNameOutCovSst2,        fileNameOutCovAcc,        fileNameOutCovPod1,        fileNameOutCovPod2;
    FileName fileNameResidualsSst1,     fileNameResidualsSst2,     fileNameResidualsAcc,     fileNameResidualsPod1,     fileNameResidualsPod2;

    FileName fileNameArcList;
    FileName fileNameInArcSigmaSst1,   fileNameInArcSigmaSst2,   fileNameInArcSigmaAcc,   fileNameInArcSigmaPod1,   fileNameInArcSigmaPod2;
    FileName fileNameInEpochSigmaSst1, fileNameInEpochSigmaSst2, fileNameInEpochSigmaAcc, fileNameInEpochSigmaPod1, fileNameInEpochSigmaPod2;
    FileName fileNameInCovSst1,        fileNameInCovSst2,        fileNameInCovAcc,        fileNameInCovPod1,        fileNameInCovPod2;

    std::vector<FileName> fileNamesCovSst1, fileNamesCovSst2, fileNamesCovAcc;
    FileName              fileNameInSigmasCovSst1, fileNameInSigmasCovSst2, fileNameInSigmasCovAcc;

    FileName    fileNameInCovEpochPod1,  fileNameInCovEpochPod2;
    Double      sigma0Sst1=1, sigma0Sst2=1, sigma0Acc=1, sigma0Pod1=1, sigma0Pod2=1;
    UInt        iterCount;
    std::string iterVariableName;
    Double      downweightPod;
    UInt        defaultBlockSize;

    AutoregressiveModelSequencePtr arSequence;
    ParameterSelectorPtr           parameterShortTime;

    readConfig(config, "outputfileSolution",      fileNameSolution, Config::OPTIONAL, "", "estimated parameter vector (static part only)");
    readConfig(config, "outputfileSigmax",        fileNameSigmax,   Config::OPTIONAL, "", "standard deviations of the parameters (sqrt of the diagonal of the inverse normal equation)");
    readConfig(config, "outputfileParameterName", fileNameParaName, Config::OPTIONAL, "", "estimated signal parameters (index is appended)");
    if(readConfigSequence(config, "estimateArcSigmas", Config::OPTIONAL, "", ""))
    {
      estimateArcSigmas = TRUE;
      readConfig(config, "outputfileSigmasPerArcSst1", fileNameOutArcSigmaSst1, Config::OPTIONAL, "", "accuracies of each arc (SST1)");
      readConfig(config, "outputfileSigmasPerArcSst2", fileNameOutArcSigmaSst2, Config::OPTIONAL, "", "accuracies of each arc (SST2)");
      readConfig(config, "outputfileSigmasPerArcAcc",  fileNameOutArcSigmaAcc,  Config::OPTIONAL, "", "accuracies of each arc (ACC)");
      readConfig(config, "outputfileSigmasPerArcPod1", fileNameOutArcSigmaPod1, Config::OPTIONAL, "", "accuracies of each arc (POD1)");
      readConfig(config, "outputfileSigmasPerArcPod2", fileNameOutArcSigmaPod2, Config::OPTIONAL, "", "accuracies of each arc (POD2)");
      endSequence(config);
    }
    if(readConfigSequence(config, "estimateEpochSigmas", Config::OPTIONAL, "", ""))
    {
      estimateEpochSigmas = TRUE;
      readConfig(config, "outputfileSigmasPerEpochSst1", fileNameOutEpochSigmaSst1, Config::OPTIONAL, "", "accuracies of each epoch (SST1)");
      readConfig(config, "outputfileSigmasPerEpochSst2", fileNameOutEpochSigmaSst2, Config::OPTIONAL, "", "accuracies of each epoch (SST2)");
      readConfig(config, "outputfileSigmasPerEpochAcc",  fileNameOutEpochSigmaAcc,  Config::OPTIONAL, "", "accuracies of each epoch (ACC)");
      readConfig(config, "outputfileSigmasPerEpochPod1", fileNameOutEpochSigmaPod1, Config::OPTIONAL, "", "accuracies of each epoch (POD1)");
      readConfig(config, "outputfileSigmasPerEpochPod2", fileNameOutEpochSigmaPod2, Config::OPTIONAL, "", "accuracies of each epoch (POD2)");
      endSequence(config);
    }
    if(readConfigSequence(config, "estimateCovarianceFunctions", Config::OPTIONAL, "", ""))
    {
      estimateCovarianceFunctionVCE = TRUE;
      readConfig(config, "outputfileCovarianceFunctionSst1", fileNameOutCovSst1, Config::OPTIONAL, "", "covariance function");
      readConfig(config, "outputfileCovarianceFunctionSst2", fileNameOutCovSst2, Config::OPTIONAL, "", "covariance function");
      readConfig(config, "outputfileCovarianceFunctionAcc",  fileNameOutCovAcc,  Config::OPTIONAL, "", "covariance function");
      readConfig(config, "outputfileCovarianceFunctionPod1", fileNameOutCovPod1, Config::OPTIONAL, "", "covariance functions for along, cross, radial direction");
      readConfig(config, "outputfileCovarianceFunctionPod2", fileNameOutCovPod2, Config::OPTIONAL, "", "covariance functions for along, cross, radial direction");
      endSequence(config);
    }
    if(readConfigSequence(config, "computeResiduals", Config::OPTIONAL, "", ""))
    {
      estimateResiduals = TRUE;
      readConfig(config, "outputfileSst1Residuals", fileNameResidualsSst1, Config::OPTIONAL, "", "");
      readConfig(config, "outputfileSst2Residuals", fileNameResidualsSst2, Config::OPTIONAL, "", "");
      readConfig(config, "outputfileAccResiduals",  fileNameResidualsAcc,  Config::OPTIONAL, "", "");
      readConfig(config, "outputfilePod1Residuals", fileNameResidualsPod1, Config::OPTIONAL, "", "");
      readConfig(config, "outputfilePod2Residuals", fileNameResidualsPod2, Config::OPTIONAL, "", "");
      endSequence(config);
    }
    readConfig(config, "observation", observationMisc, Config::MUSTSET,  "", "");
    if(readConfigSequence(config, "covarianceSst1", Config::MUSTSET, "", ""))
    {
      readConfig(config, "sigma",                              sigma0Sst1,               Config::DEFAULT,  "1", "apriori factor of covariance function");
      readConfig(config, "inputfileSigmasPerArc",              fileNameInArcSigmaSst1,   Config::OPTIONAL, "",  "apriori different accuaries for each arc (multiplicated with sigma)");
      readConfig(config, "inputfileSigmasPerEpoch",            fileNameInEpochSigmaSst1, Config::OPTIONAL, "",  "apriori different accuaries for each epoch");
      readConfig(config, "inputfileCovarianceFunction",        fileNameInCovSst1,        Config::OPTIONAL, "",  "approximate covariances in time");
      readConfig(config, "inputfileCovarianceMatrixArc",       fileNamesCovSst1,         Config::OPTIONAL, "",  "Must be given per sst arc with correct dimensions.");
      readConfig(config, "inputfileSigmasCovarianceMatrixArc", fileNameInSigmasCovSst1,  Config::OPTIONAL, "",  "Vector with one sigma for each <inputfileCovarianceMatrixArc>");
      readConfig(config, "sampling",                           samplingSst,              Config::DEFAULT,  "5", "[seconds] sampling of the covariance function");
      endSequence(config);
    }
    if(readConfigSequence(config, "covarianceSst2", Config::MUSTSET, "", ""))
    {
      readConfig(config, "sigma",                              sigma0Sst2,               Config::DEFAULT,  "1", "apriori factor of covariance function");
      readConfig(config, "inputfileSigmasPerArc",              fileNameInArcSigmaSst2,   Config::OPTIONAL, "",  "apriori different accuaries for each arc (multiplicated with sigma)");
      readConfig(config, "inputfileSigmasPerEpoch",            fileNameInEpochSigmaSst2, Config::OPTIONAL, "",  "apriori different accuaries for each epoch");
      readConfig(config, "inputfileCovarianceFunction",        fileNameInCovSst2,        Config::OPTIONAL, "",  "approximate covariances in time");
      readConfig(config, "inputfileCovarianceMatrixArc",       fileNamesCovSst2,         Config::OPTIONAL, "",  "Must be given per sst arc with correct dimensions.");
      readConfig(config, "inputfileSigmasCovarianceMatrixArc", fileNameInSigmasCovSst2,  Config::OPTIONAL, "",  "Vector with one sigma for each <inputfileCovarianceMatrixArc>");
      endSequence(config);
    }
    if(readConfigSequence(config, "covarianceAcc", Config::MUSTSET, "", ""))
    {
      readConfig(config, "sigma",                              sigma0Acc,                Config::DEFAULT,  "1", "apriori factor of covariance function");
      readConfig(config, "inputfileSigmasPerArc",              fileNameInArcSigmaAcc,    Config::OPTIONAL, "",  "apriori different accuaries for each arc (multplicated with sigma)");
      readConfig(config, "inputfileSigmasPerEpoch",            fileNameInEpochSigmaAcc,  Config::OPTIONAL, "",  "apriori different accuaries for each epoch");
      readConfig(config, "inputfileCovarianceFunction",        fileNameInCovAcc,         Config::OPTIONAL, "",  "approximate covariances in time");
      readConfig(config, "inputfileCovarianceMatrixArc",       fileNamesCovAcc,          Config::OPTIONAL, "",  "Must be given per sst arc with correct dimensions.");
      readConfig(config, "inputfileSigmasCovarianceMatrixArc", fileNameInSigmasCovAcc,   Config::OPTIONAL, "",  "Vector with one sigma for each <inputfileCovarianceMatrixArc>");
      endSequence(config);
    }
    if(readConfigSequence(config, "covariancePod1", Config::MUSTSET, "", ""))
    {
      readConfig(config, "sigma",                        sigma0Pod1,               Config::DEFAULT,  "1",  "apriori factor of covariance function");
      readConfig(config, "inputfileSigmasPerArc",        fileNameInArcSigmaPod1,   Config::OPTIONAL, "",   "apriori different accuaries for each arc (multiplicated with sigma)");
      readConfig(config, "inputfileSigmasPerEpoch",      fileNameInEpochSigmaPod1, Config::OPTIONAL, "",  "apriori different accuaries for each epoch");
      readConfig(config, "inputfileCovarianceFunction",  fileNameInCovPod1,        Config::OPTIONAL, "",   "approximate covariances in time");
      readConfig(config, "inputfileCovariancePodEpoch",  fileNameInCovEpochPod1,   Config::OPTIONAL, "",   "3x3 epoch covariances");
      readConfig(config, "sampling",                     samplingPod1,             Config::DEFAULT,  "30", "[seconds] sampling of the covariance function");
      endSequence(config);
    }
    if(readConfigSequence(config, "covariancePod2", Config::MUSTSET, "", ""))
    {
      readConfig(config, "sigma",                        sigma0Pod2,               Config::DEFAULT,  "1",  "apriori factor of covariance function");
      readConfig(config, "inputfileSigmasPerArc",        fileNameInArcSigmaPod2,   Config::OPTIONAL, "",   "apriori different accuaries for each arc (multiplicated with sigma)");
      readConfig(config, "inputfileSigmasPerEpoch",      fileNameInEpochSigmaPod2, Config::OPTIONAL, "",   "apriori different accuaries for each epoch");
      readConfig(config, "inputfileCovarianceFunction",  fileNameInCovPod2,        Config::OPTIONAL, "",   "approximate covariances in time");
      readConfig(config, "inputfileCovariancePodEpoch",  fileNameInCovEpochPod2,   Config::OPTIONAL, "",   "3x3 epoch covariances");
      readConfig(config, "sampling",                     samplingPod2,             Config::DEFAULT,  "30", "[seconds] sampling of the covariance function");
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
    sigmaSst1 = sigmaSst2 = sigmaAcc = sigmaPod1 = sigmaPod2 = Vector(arcCount);
    for(UInt arcNo=0; arcNo<arcCount; arcNo++)
      sigmaSst1(arcNo) = sigmaSst2(arcNo) = sigmaAcc(arcNo) = sigmaPod1(arcNo) = sigmaPod2(arcNo) = 1.0;
    if(!fileNameInArcSigmaSst1.empty()) readFileMatrix(fileNameInArcSigmaSst1, sigmaSst1);
    if(!fileNameInArcSigmaSst2.empty()) readFileMatrix(fileNameInArcSigmaSst2, sigmaSst2);
    if(!fileNameInArcSigmaSst2.empty()) readFileMatrix(fileNameInArcSigmaAcc,  sigmaAcc);
    if(!fileNameInArcSigmaPod1.empty()) readFileMatrix(fileNameInArcSigmaPod1, sigmaPod1);
    if(!fileNameInArcSigmaPod1.empty()) readFileMatrix(fileNameInArcSigmaPod2, sigmaPod2);
    if(sigmaSst1.rows() != arcCount) throw(Exception("sigmasPerArc (SST1) contains wrong number of arcs"));
    if(sigmaSst2.rows() != arcCount) throw(Exception("sigmasPerArc (SST2) contains wrong number of arcs"));
    if(sigmaAcc.rows()  != arcCount) throw(Exception("sigmasPerArc (ACC)  contains wrong number of arcs"));
    if(sigmaPod1.rows() != arcCount) throw(Exception("sigmasPerArc (POD1) contains wrong number of arcs"));
    if(sigmaPod2.rows() != arcCount) throw(Exception("sigmasPerArc (POD2) contains wrong number of arcs"));

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
    CovSst1.resize(arcCount);
    CovSst2.resize(arcCount);
    CovAcc.resize(arcCount);
    if(fileNamesCovSst1.size() || fileNamesCovSst2.size() || fileNamesCovAcc.size())
    {
      logStatus<<"read arc-wise sst covariance matrices"<<Log::endl;
      Parallel::forEachProcess(arcCount, [&](UInt arcNo)
      {
        CovSst1.at(arcNo).resize(fileNamesCovSst1.size());
        for(UInt i=0; i<fileNamesCovSst1.size(); i++)
          readFileMatrix(fileNamesCovSst1.at(i).appendBaseName(".arc"+arcNo%"%03i"s), CovSst1.at(arcNo).at(i));

        CovSst2.at(arcNo).resize(fileNamesCovSst2.size());
        for(UInt i=0; i<fileNamesCovSst2.size(); i++)
          readFileMatrix(fileNamesCovSst2.at(i).appendBaseName(".arc"+arcNo%"%03i"s), CovSst2.at(arcNo).at(i));

        CovAcc.at(arcNo).resize (fileNamesCovAcc.size());
        for(UInt i=0; i<fileNamesCovAcc.size(); i++)
          readFileMatrix(fileNamesCovAcc.at(i).appendBaseName(".arc"+arcNo%"%03i"s),  CovAcc.at(arcNo).at(i));
      }, processNo, comm);
    }

    sigmaCovSst1 = Vector(fileNamesCovSst1.size(), 1.);
    sigmaCovSst2 = Vector(fileNamesCovSst2.size(), 1.);
    sigmaCovAcc  = Vector(fileNamesCovAcc.size(),  1.);
    if(!fileNameInSigmasCovSst1.empty()) readFileMatrix(fileNameInSigmasCovSst1, sigmaCovSst1);
    if(!fileNameInSigmasCovSst2.empty()) readFileMatrix(fileNameInSigmasCovSst2, sigmaCovSst2);
    if(!fileNameInSigmasCovAcc.empty())  readFileMatrix(fileNameInSigmasCovAcc,  sigmaCovAcc);

    if((sigmaCovSst1.rows() != fileNamesCovSst1.size()) || (sigmaCovSst2.rows() != fileNamesCovSst2.size()) || (sigmaCovAcc.rows() != fileNamesCovAcc.size()))
      throw(Exception("Number of sigmas not compatible with number of given arc-wise SST covariance matrices"));

    // =============================================

    // Determine max. length of ovariance functions
    // --------------------------------------------
    UInt covLengthSst  = 0;
    UInt covLengthPod1 = 0;
    UInt covLengthPod2 = 0;
    for(UInt arcNo=0; arcNo<arcCount; arcNo++)
    {
      covLengthSst  = std::max(covLengthSst,  observationArc.at(arcNo).timesSst.size());
      covLengthPod1 = std::max(covLengthPod1, observationArc.at(arcNo).timesPod1.size() ? static_cast<UInt>(std::round((observationArc.at(arcNo).timesPod1.back() - observationArc.at(arcNo).timesPod1.front()).seconds()/samplingPod1) + 1) : 0);
      covLengthPod2 = std::max(covLengthPod2, observationArc.at(arcNo).timesPod2.size() ? static_cast<UInt>(std::round((observationArc.at(arcNo).timesPod2.back() - observationArc.at(arcNo).timesPod2.front()).seconds()/samplingPod2) + 1) : 0);
    }
    Parallel::reduceMax(covLengthSst,  0, comm);  Parallel::broadCast(covLengthSst, 0, comm);
    Parallel::reduceMax(covLengthPod1, 0, comm); Parallel::broadCast(covLengthPod1, 0, comm);
    Parallel::reduceMax(covLengthPod2, 0, comm); Parallel::broadCast(covLengthPod2, 0, comm);

    // ===================================================

    // transformation PSD <-> covFunc
    // ------------------------------
    CosTransformSst  = Vce::cosTransform(covLengthSst);
    CosTransformPod1 = Vce::cosTransform(covLengthPod1);
    CosTransformPod2 = Vce::cosTransform(covLengthPod2);

    // init covariance function
    // ------------------------
    covFuncSst1 = Vce::readCovarianceFunction(fileNameInCovSst1, covLengthSst,  1, samplingSst);
    covFuncSst2 = Vce::readCovarianceFunction(fileNameInCovSst2, covLengthSst,  1, samplingSst);
    covFuncAcc  = Vce::readCovarianceFunction(fileNameInCovAcc,  covLengthSst,  1, samplingSst);
    covFuncPod1 = Vce::readCovarianceFunction(fileNameInCovPod1, covLengthPod1, 3, samplingPod1);
    covFuncPod2 = Vce::readCovarianceFunction(fileNameInCovPod2, covLengthPod2, 3, samplingPod2);
    covFuncSst1.column(1,1) *= std::pow(sigma0Sst1, 2);
    covFuncSst2.column(1,1) *= std::pow(sigma0Sst2, 2);
    covFuncAcc.column(1,1)  *= std::pow(sigma0Acc,  2);
    covFuncPod1.column(1,3) *= std::pow(sigma0Pod1, 2);
    covFuncPod2.column(1,3) *= std::pow(sigma0Pod2, 2);
    PsdSst1 = CosTransformSst  * covFuncSst1.column(1,1);
    PsdSst2 = CosTransformSst  * covFuncSst2.column(1,1);
    PsdAcc  = CosTransformSst  * covFuncAcc.column(1,1);
    PsdPod1 = CosTransformPod1 * covFuncPod1.column(1,3);
    PsdPod2 = CosTransformPod2 * covFuncPod2.column(1,3);

    // =============================================

    // init epoch sigmas
    // -----------------
    arcListEpochSigmaSst1.resize(arcCount);
    arcListEpochSigmaSst2.resize(arcCount);
    arcListEpochSigmaAcc.resize(arcCount);
    arcListEpochSigmaPod1.resize(arcCount);
    arcListEpochSigmaPod2.resize(arcCount);

    if(estimateEpochSigmas)
    {
      InstrumentFile fileSst1(fileNameInEpochSigmaSst1);
      InstrumentFile fileSst2(fileNameInEpochSigmaSst2);
      InstrumentFile fileAcc(fileNameInEpochSigmaAcc);
      InstrumentFile filePod1(fileNameInEpochSigmaPod1);
      InstrumentFile filePod2(fileNameInEpochSigmaPod2);

      for(UInt arcNo=0; arcNo<arcCount; arcNo++)
        if(Parallel::myRank(comm) == processNo.at(arcNo))
        {
          arcListEpochSigmaSst1.at(arcNo) = fileSst1.readArc(arcNo);
          arcListEpochSigmaSst2.at(arcNo) = fileSst2.readArc(arcNo);
          arcListEpochSigmaAcc.at(arcNo)  = fileAcc.readArc(arcNo);
          arcListEpochSigmaPod1.at(arcNo) = filePod1.readArc(arcNo);
          arcListEpochSigmaPod2.at(arcNo) = filePod2.readArc(arcNo);

          if(arcListEpochSigmaSst1.at(arcNo).size()==0)
            for(UInt i=0; i<observationArc.at(arcNo).timesSst.size(); i++)
            {
              ObservationSigmaEpoch epoch;
              epoch.time = observationArc.at(arcNo).timesSst.at(i);
              arcListEpochSigmaSst1.at(arcNo).push_back(epoch);
            }

          if(arcListEpochSigmaSst2.at(arcNo).size()==0)
            for(UInt i=0; i<observationArc.at(arcNo).timesSst.size(); i++)
            {
              ObservationSigmaEpoch epoch;
              epoch.time = observationArc.at(arcNo).timesSst.at(i);
              arcListEpochSigmaSst2.at(arcNo).push_back(epoch);
            }

          if(arcListEpochSigmaAcc.at(arcNo).size()==0)
            for(UInt i=0; i<observationArc.at(arcNo).timesSst.size(); i++)
            {
              ObservationSigmaEpoch epoch;
              epoch.time = observationArc.at(arcNo).timesSst.at(i);
              arcListEpochSigmaAcc.at(arcNo).push_back(epoch);
            }

          if(arcListEpochSigmaPod1.at(arcNo).size()==0)
            for(UInt i=0; i<observationArc.at(arcNo).timesPod1.size(); i++)
            {
              ObservationSigmaEpoch epoch;
              epoch.time = observationArc.at(arcNo).timesPod1.at(i);
              arcListEpochSigmaPod1.at(arcNo).push_back(epoch);
            }

          if(arcListEpochSigmaPod2.at(arcNo).size()==0)
            for(UInt i=0; i<observationArc.at(arcNo).timesPod2.size(); i++)
            {
              ObservationSigmaEpoch epoch;
              epoch.time = observationArc.at(arcNo).timesPod2.at(i);
              arcListEpochSigmaPod2.at(arcNo).push_back(epoch);
            }
        }
    } // if(estimateEpochSigmas)

    // =============================================

    // Iteration
    // ---------
    VariableList variableIteration;
    if(iterVariableName.size())
      addVariable(iterVariableName, variableIteration);

    Double sigma2ShortTimeModel = 1.;
    for(UInt iter=0; iter<iterCount; iter++)
    {
      logInfo<<"starting iteration "<<iter<<Log::endl;
      if(iterVariableName.size())
        variableIteration[iterVariableName]->setValue(iter);

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
        Double sigma = normals.solve(x, Wz);
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
          if((s2==s2) && (s2>0))
            sigma2ShortTimeModel = s2;
        }
      } // if(countAParameter)
      Parallel::barrier(comm);

      // =============================================

      if(estimateResiduals || estimateEpochSigmas)
      {
        logStatus<<"compute residuals"<<Log::endl;
        arcListResidualsSst1.clear(); arcListResidualsSst1.resize(arcCount);
        arcListResidualsSst2.clear(); arcListResidualsSst2.resize(arcCount);
        arcListResidualsAcc.clear();  arcListResidualsAcc.resize(arcCount);
        arcListResidualsPod1.clear(); arcListResidualsPod1.resize(arcCount);
        arcListResidualsPod2.clear(); arcListResidualsPod2.resize(arcCount);
        Parallel::forEachProcess(arcCount, [this](UInt arcNo) {computeResiduals(arcNo);}, processNo, comm);
        Parallel::forEachProcess(arcListResidualsSst1, [this](UInt arcNo) {return arcListResidualsSst1.at(arcNo);}, processNo, comm);
        Parallel::forEachProcess(arcListResidualsSst2, [this](UInt arcNo) {return arcListResidualsSst2.at(arcNo);}, processNo, comm);
        Parallel::forEachProcess(arcListResidualsAcc,  [this](UInt arcNo) {return arcListResidualsAcc.at(arcNo);},  processNo, comm);
        Parallel::forEachProcess(arcListResidualsPod1, [this](UInt arcNo) {return arcListResidualsPod1.at(arcNo);}, processNo, comm);
        Parallel::forEachProcess(arcListResidualsPod2, [this](UInt arcNo) {return arcListResidualsPod2.at(arcNo);}, processNo, comm);

        if(Parallel::isMaster(comm) && (!fileNameResidualsSst1.empty()))
        {
          logStatus<<"write residual (SST1) file <"<<fileNameResidualsSst1(variableIteration)<<">"<<Log::endl;
          InstrumentFile::write(fileNameResidualsSst1(variableIteration), arcListResidualsSst1);
        }

        if(Parallel::isMaster(comm) && (!fileNameResidualsSst2.empty()))
        {
          logStatus<<"write residual (SST2) file <"<<fileNameResidualsSst2(variableIteration)<<">"<<Log::endl;
          InstrumentFile::write(fileNameResidualsSst2(variableIteration), arcListResidualsSst2);
        }

        if(Parallel::isMaster(comm) && (!fileNameResidualsAcc.empty()))
        {
          logStatus<<"write residual (ACC) file <"<<fileNameResidualsAcc(variableIteration)<<">"<<Log::endl;
          InstrumentFile::write(fileNameResidualsAcc(variableIteration), arcListResidualsAcc);
        }

        if(Parallel::isMaster(comm) && (!fileNameResidualsPod1.empty()))
        {
          logStatus<<"write residual (POD1) file <"<<fileNameResidualsPod1(variableIteration)<<">"<<Log::endl;
          InstrumentFile::write(fileNameResidualsPod1(variableIteration), arcListResidualsPod1);
        }

        if(Parallel::isMaster(comm) && (!fileNameResidualsPod2.empty()))
        {
          logStatus<<"write residual (POD2) file <"<<fileNameResidualsPod2(variableIteration)<<">"<<Log::endl;
          InstrumentFile::write(fileNameResidualsPod2(variableIteration), arcListResidualsPod2);
        }
      } // if(estimateResiduals)

      // =============================================

      // compute redundancies
      // --------------------
      if(estimateArcSigmas || estimateCovarianceFunctionVCE)
      {
        logStatus<<"compute redundancies"<<Log::endl;
        sigmaSst1New = sigmaSst2New = sigmaAccNew = sigmaPod1New = sigmaPod2New = Vector(arcCount);
        ePeSst1 = redundancySst1 = Matrix(covLengthSst,  1);
        ePeSst2 = redundancySst2 = Matrix(covLengthSst,  1);
        ePeAcc  = redundancyAcc  = Matrix(covLengthSst,  1);
        ePePod1 = redundancyPod1 = Matrix(covLengthPod1, 3); // for x,y,z
        ePePod2 = redundancyPod2 = Matrix(covLengthPod2, 3); // for x,y,z
        Parallel::forEachProcess(arcCount, [this](UInt arcNo) {computeRedundancies(arcNo);}, processNo, comm);
      }

      // =============================================

      // sigmas per arc
      // --------------
      if(estimateArcSigmas)
      {
        Parallel::reduceSum(sigmaSst1New, 0, comm);
        Parallel::reduceSum(sigmaSst2New, 0, comm);
        Parallel::reduceSum(sigmaAccNew,  0, comm);
        Parallel::reduceSum(sigmaPod1New, 0, comm);
        Parallel::reduceSum(sigmaPod2New, 0, comm);
        if(Parallel::isMaster(comm))
        {
          sigmaSst1 = sigmaSst1New;
          sigmaSst2 = sigmaSst2New;
          sigmaAcc  = sigmaAccNew;
          sigmaPod1 = sigmaPod1New;
          sigmaPod2 = sigmaPod2New;

          Double sigma0Sst1 = Vce::meanSigma(sigmaSst1);
          Double sigma0Sst2 = Vce::meanSigma(sigmaSst2);
          Double sigma0Acc  = Vce::meanSigma(sigmaAcc);
          Double sigma0Pod1 = Vce::meanSigma(sigmaPod1);
          Double sigma0Pod2 = Vce::meanSigma(sigmaPod2);

          sigmaSst1 *= 1./sigma0Sst1;
          sigmaSst2 *= 1./sigma0Sst2;
          sigmaAcc  *= 1./sigma0Acc;
          sigmaPod1 *= 1./sigma0Pod1;
          sigmaPod2 *= 1./sigma0Pod2;

          logInfo<<"  SST1 sigma per arc (median): "<<sigma0Sst1<<Log::endl;
          logInfo<<"  SST2 sigma per arc (median): "<<sigma0Sst2<<Log::endl;
          logInfo<<"  ACC  sigma per arc (median): "<<sigma0Acc<<Log::endl;
          logInfo<<"  POD1 sigma per arc (median): "<<sigma0Pod1<<Log::endl;
          logInfo<<"  POD2 sigma per arc (median): "<<sigma0Pod2<<Log::endl;

          if(!fileNameOutArcSigmaSst1.empty())
          {
            logStatus<<"write arc sigma file <"<<fileNameOutArcSigmaSst1(variableIteration)<<">"<<Log::endl;
            writeFileMatrix(fileNameOutArcSigmaSst1(variableIteration), sigmaSst1);
          }
          if(!fileNameOutArcSigmaSst2.empty())
          {
            logStatus<<"write arc sigma file <"<<fileNameOutArcSigmaSst2(variableIteration)<<">"<<Log::endl;
            writeFileMatrix(fileNameOutArcSigmaSst2(variableIteration), sigmaSst2);
          }
          if(!fileNameOutArcSigmaAcc.empty())
          {
            logStatus<<"write arc sigma file <"<<fileNameOutArcSigmaAcc(variableIteration)<<">"<<Log::endl;
            writeFileMatrix(fileNameOutArcSigmaAcc(variableIteration), sigmaAcc);
          }
          if(!fileNameOutArcSigmaPod1.empty())
          {
            logStatus<<"write arc sigma file <"<<fileNameOutArcSigmaPod1(variableIteration)<<">"<<Log::endl;
            writeFileMatrix(fileNameOutArcSigmaPod1(variableIteration), sigmaPod1);
          }
          if(!fileNameOutArcSigmaPod2.empty())
          {
            logStatus<<"write arc sigma file <"<<fileNameOutArcSigmaPod2(variableIteration)<<">"<<Log::endl;
            writeFileMatrix(fileNameOutArcSigmaPod2(variableIteration), sigmaPod2);
          }
        }
        Parallel::broadCast(sigmaSst1, 0, comm);
        Parallel::broadCast(sigmaSst2, 0, comm);
        Parallel::broadCast(sigmaAcc,  0, comm);
        Parallel::broadCast(sigmaPod1, 0, comm);
        Parallel::broadCast(sigmaPod2, 0, comm);
      } // if(estimateArcSigmas)

      // =============================================

      if(estimateEpochSigmas)
      {
        logStatus<<"compute epoch sigmas"<<Log::endl;
        Parallel::forEachProcess(arcCount, [this](UInt arcNo) {computeEpochSigmas(arcNo);}, processNo, comm);

        logStatus<<"collect epoch sigmas"<<Log::endl;
        Parallel::forEachProcess(arcListEpochSigmaSst1, [this](UInt arcNo) {return arcListEpochSigmaSst1.at(arcNo);}, processNo, comm);
        Parallel::forEachProcess(arcListEpochSigmaSst2, [this](UInt arcNo) {return arcListEpochSigmaSst2.at(arcNo);}, processNo, comm);
        Parallel::forEachProcess(arcListEpochSigmaAcc,  [this](UInt arcNo) {return arcListEpochSigmaAcc.at(arcNo);},  processNo, comm);
        Parallel::forEachProcess(arcListEpochSigmaPod1, [this](UInt arcNo) {return arcListEpochSigmaPod1.at(arcNo);}, processNo, comm);
        Parallel::forEachProcess(arcListEpochSigmaPod2, [this](UInt arcNo) {return arcListEpochSigmaPod2.at(arcNo);}, processNo, comm);

        if(Parallel::isMaster(comm) && (!fileNameOutEpochSigmaSst1.empty()))
        {
          logStatus<<"write epoch sigma (SST1) file <"<<fileNameOutEpochSigmaSst1(variableIteration)<<">"<<Log::endl;
          InstrumentFile::write(fileNameOutEpochSigmaSst1(variableIteration), arcListEpochSigmaSst1);
        }

        if(Parallel::isMaster(comm) && (!fileNameOutEpochSigmaSst2.empty()))
        {
          logStatus<<"write epoch sigma (SST2) file <"<<fileNameOutEpochSigmaSst2(variableIteration)<<">"<<Log::endl;
          InstrumentFile::write(fileNameOutEpochSigmaSst2(variableIteration), arcListEpochSigmaSst2);
        }

        if(Parallel::isMaster(comm) && (!fileNameOutEpochSigmaAcc.empty()))
        {
          logStatus<<"write epoch sigma (Acc) file <"<<fileNameOutEpochSigmaAcc(variableIteration)<<">"<<Log::endl;
          InstrumentFile::write(fileNameOutEpochSigmaAcc(variableIteration), arcListEpochSigmaAcc);
        }

        if(Parallel::isMaster(comm) && (!fileNameOutEpochSigmaPod1.empty()))
        {
          logStatus<<"write epoch sigma (POD1) file <"<<fileNameOutEpochSigmaPod1(variableIteration)<<">"<<Log::endl;
          InstrumentFile::write(fileNameOutEpochSigmaPod1(variableIteration), arcListEpochSigmaPod1);
        }

        if(Parallel::isMaster(comm) && (!fileNameOutEpochSigmaPod2.empty()))
        {
          logStatus<<"write epoch sigma (POD2) file <"<<fileNameOutEpochSigmaPod2(variableIteration)<<">"<<Log::endl;
          InstrumentFile::write(fileNameOutEpochSigmaPod2(variableIteration), arcListEpochSigmaPod2);
        }

        Parallel::barrier(comm);
      } // if(estimateEpochSigmas)

      // =============================================

      // estimate new PSD through variance component estimation
      // ------------------------------------------------------
      if(estimateCovarianceFunctionVCE)
      {
        Parallel::reduceSum(ePeSst1, 0, comm);
        Parallel::reduceSum(ePeSst2, 0, comm);
        Parallel::reduceSum(ePeAcc,  0, comm);
        Parallel::reduceSum(ePePod1, 0, comm);
        Parallel::reduceSum(ePePod2, 0, comm);
        Parallel::reduceSum(redundancySst1, 0, comm);
        Parallel::reduceSum(redundancySst2, 0, comm);
        Parallel::reduceSum(redundancyAcc,  0, comm);
        Parallel::reduceSum(redundancyPod1, 0, comm);
        Parallel::reduceSum(redundancyPod2, 0, comm);

        logStatus<<"compute psd through variance component estimation"<<Log::endl;
        if(Parallel::isMaster(comm))
        {
          Double maxFactorSst1 = 0;
          Double maxFactorSst2 = 0;
          Double maxFactorAcc  = 0;
          Double maxFactorPod  = 0;
          Vce::estimatePsd(ePeSst1, redundancySst1, PsdSst1, maxFactorSst1);
          Vce::estimatePsd(ePeSst2, redundancySst2, PsdSst2, maxFactorSst2);
          Vce::estimatePsd(ePeAcc,  redundancyAcc,  PsdAcc,  maxFactorAcc);
          Vce::estimatePsd(ePePod1, redundancyPod1, PsdPod1, maxFactorPod);
          Vce::estimatePsd(ePePod2, redundancyPod2, PsdPod2, maxFactorPod);

          maxFactorPod /= downweightPod;
          logInfo<<"  max. PSD adjustment factor (SST1): "<<maxFactorSst1<<Log::endl;
          logInfo<<"  max. PSD adjustment factor (SST2): "<<maxFactorSst2<<Log::endl;
          logInfo<<"  max. PSD adjustment factor (ACC):  "<<maxFactorAcc<<Log::endl;
          logInfo<<"  max. PSD adjustment factor (POD):  "<<maxFactorPod<<Log::endl;

        } // if(Parallel::isMaster(comm))
        Parallel::broadCast(PsdSst1, 0, comm);
        Parallel::broadCast(PsdSst2, 0, comm);
        Parallel::broadCast(PsdAcc,  0, comm);
        Parallel::broadCast(PsdPod1, 0, comm);
        Parallel::broadCast(PsdPod2, 0, comm);
        // compute new covariance function
        copy(CosTransformSst  * PsdSst1, covFuncSst1.column(1,1));
        copy(CosTransformSst  * PsdSst2, covFuncSst2.column(1,1));
        copy(CosTransformSst  * PsdAcc,  covFuncAcc.column(1,1));
        copy(CosTransformPod1 * PsdPod1, covFuncPod1.column(1,3));
        copy(CosTransformPod2 * PsdPod2, covFuncPod2.column(1,3));
      } // if(estimateCovarianceFunctions)

      // =============================================

      // Write covariance function to file
      // ---------------------------------
      if(estimateCovarianceFunctionVCE)
      {
        if(Parallel::isMaster(comm) && !fileNameOutCovSst1.empty())
        {
          logStatus<<"write covariance function file <"<<fileNameOutCovSst1(variableIteration)<<">"<<Log::endl;
          writeFileMatrix(fileNameOutCovSst1(variableIteration), covFuncSst1);
        }
        if(Parallel::isMaster(comm) && !fileNameOutCovSst2.empty())
        {
          logStatus<<"write covariance function file <"<<fileNameOutCovSst2(variableIteration)<<">"<<Log::endl;
          writeFileMatrix(fileNameOutCovSst2(variableIteration), covFuncSst2);
        }
        if(Parallel::isMaster(comm) && !fileNameOutCovAcc.empty())
        {
          logStatus<<"write covariance function file <"<<fileNameOutCovAcc(variableIteration)<<">"<<Log::endl;
          writeFileMatrix(fileNameOutCovAcc(variableIteration), covFuncAcc);
        }
        if(Parallel::isMaster(comm) && !fileNameOutCovPod1.empty())
        {
          logStatus<<"write covariance function file <"<<fileNameOutCovPod1(variableIteration)<<">"<<Log::endl;
          writeFileMatrix(fileNameOutCovPod1(variableIteration), covFuncPod1);
        }
        if(Parallel::isMaster(comm) && !fileNameOutCovPod2.empty())
        {
          logStatus<<"write covariance function file <"<<fileNameOutCovPod2(variableIteration)<<">"<<Log::endl;
          writeFileMatrix(fileNameOutCovPod2(variableIteration), covFuncPod2);
        }
      }

      PsdPod1                 *= std::pow(downweightPod,2);
      PsdPod2                 *= std::pow(downweightPod,2);
      covFuncPod1.column(1,3) *= std::pow(downweightPod,2);
      covFuncPod2.column(1,3) *= std::pow(downweightPod,2);
    } // for(iter)

    // =============================================
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

void PreprocessingDualSst::decorrelate(UInt arcNo, UInt countSst, UInt countPod1, UInt countPod2,
                                       Matrix &CovSst1, Matrix &CovSst2, Matrix &CovAcc,
                                       Matrix &WSst, Matrix &WPod1, Matrix &WPod2, const std::list<MatrixSlice> &A)
{
  try
  {
    UInt obsCount = 0;
    const UInt idxSst  = obsCount; obsCount += 2*countSst;
    const UInt idxPod1 = obsCount; obsCount += 3*countPod1;
    const UInt idxPod2 = obsCount; obsCount += 3*countPod2;

    if(countSst)
    {
      std::list<MatrixSlice> WA;
      for(const MatrixSlice &a : A)
        WA.push_back(a.row(idxSst, 2*countSst));

      if(this->CovSst1.at(arcNo).size())
      {
        CovSst1 = std::pow(sigmaCovSst1.at(0), 2) * this->CovSst1.at(arcNo).at(0);
        for(UInt i=1; i<this->CovSst1.at(arcNo).size(); i++)
          axpy(std::pow(sigmaCovSst1.at(i), 2), this->CovSst1.at(arcNo).at(i), CovSst1);
      }

      if(this->CovSst2.at(arcNo).size())
      {
        CovSst2 = std::pow(sigmaCovSst2.at(0), 2) * this->CovSst2.at(arcNo).at(0);
        for(UInt i=1; i<this->CovSst2.at(arcNo).size(); i++)
          axpy(std::pow(sigmaCovSst2.at(i), 2), this->CovSst2.at(arcNo).at(i), CovSst2);
      }

      if(this->CovAcc.at(arcNo).size())
      {
        CovAcc = std::pow(sigmaCovAcc.at(0), 2) * this->CovAcc.at(arcNo).at(0);
        for(UInt i=1; i<this->CovAcc.at(arcNo).size(); i++)
          axpy(std::pow(sigmaCovAcc.at(i), 2), this->CovAcc.at(arcNo).at(i), CovAcc);
      }

      CovarianceSst::covariance(observationArc.at(arcNo).timesSst, sigmaSst1(arcNo), arcListEpochSigmaSst1.at(arcNo), covFuncSst1, CovSst1);
      CovarianceSst::covariance(observationArc.at(arcNo).timesSst, sigmaSst2(arcNo), arcListEpochSigmaSst2.at(arcNo), covFuncSst2, CovSst2);
      CovarianceSst::covariance(observationArc.at(arcNo).timesSst, sigmaAcc(arcNo),  arcListEpochSigmaAcc.at(arcNo),  covFuncAcc,  CovAcc);
      WSst = observationMisc->decorrelate(CovSst1, CovSst2, CovAcc, WA);
    }

    if(countPod1)
    {
      std::list<MatrixSlice> WA;
      for(const MatrixSlice &a : A)
        WA.push_back(a.row(idxPod1, 3*countPod1));
      WPod1 = CovariancePod::decorrelate(observationArc.at(arcNo).pod1, sigmaPod1(arcNo), arcListEpochSigmaPod1.at(arcNo), fileCovEpochPod1.readArc(arcNo), covFuncPod1, WA);
    }

    if(countPod2)
    {
      std::list<MatrixSlice> WA;
      for(const MatrixSlice &a : A)
        WA.push_back(a.row(idxPod2, 3*countPod2));
      WPod2 = CovariancePod::decorrelate(observationArc.at(arcNo).pod2, sigmaPod2(arcNo), arcListEpochSigmaPod2.at(arcNo), fileCovEpochPod2.readArc(arcNo), covFuncPod2, WA);
    }
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
    Matrix CovSst1, CovSst2, CovAcc, WSst, WPod1, WPod2;
    decorrelate(arcNo, observationArc.at(arcNo).timesSst.size(),
                observationArc.at(arcNo).timesPod1.size(), observationArc.at(arcNo).timesPod2.size(),
                CovSst1, CovSst2, CovAcc, WSst, WPod1, WPod2, {Wl, WA, WB});

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
    const UInt countSst  = observationArc.at(arcNo).timesSst.size();
    const UInt countPod1 = observationArc.at(arcNo).timesPod1.size();
    const UInt countPod2 = observationArc.at(arcNo).timesPod2.size();

    UInt obsCount = 0;
    const UInt idxSst  = obsCount; obsCount += 2*countSst;
    const UInt idxPod1 = obsCount; obsCount += 3*countPod1;
    const UInt idxPod2 = obsCount; obsCount += 3*countPod2;

    // Decorrelation
    // -------------
    Matrix Wl = observationArc.at(arcNo).l;
    Matrix WA = observationArc.at(arcNo).A;
    Matrix WB = observationArc.at(arcNo).B;
    Matrix CovSst1, CovSst2, CovAcc, WSst, WPod1, WPod2;
    decorrelate(arcNo, observationArc.at(arcNo).timesSst.size(),
                observationArc.at(arcNo).timesPod1.size(), observationArc.at(arcNo).timesPod2.size(),
                CovSst1, CovSst2, CovAcc, WSst, WPod1, WPod2, {Wl, WA, WB});

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

    // variance component estimation (Sst)
    // -----------------------------------
    if(countSst)
    {
      Matrix R;
      Vector WWe;
      Vce::redundancy(WSst, We.row(idxSst,  2*countSst), WAz.row(idxSst, 2*countSst), WB.row(idxSst, 2*countSst), R, WWe);

      Double ePeSst1Sum = 0, ePeSst2Sum = 0, ePeAccSum = 0;
      Double redundancySst1Sum = 0, redundancySst2Sum = 0, redundancyAccSum = 0;

      varianceComponentEstimationPsd(R, WWe, CosTransformSst,
                                     sigmaSst1(arcNo), sigmaSst2(arcNo), sigmaAcc(arcNo),
                                     PsdSst1, PsdSst2, PsdAcc,
                                     ePeSst1, redundancySst1, ePeSst1Sum, redundancySst1Sum,
                                     ePeSst2, redundancySst2, ePeSst2Sum, redundancySst2Sum,
                                     ePeAcc,  redundancyAcc,  ePeAccSum,  redundancyAccSum);

      sigmaSst1New(arcNo) = std::sqrt(ePeSst1Sum/redundancySst1Sum) * sigmaSst1(arcNo);  // compute new sigma (for this arc)
      sigmaSst2New(arcNo) = std::sqrt(ePeSst2Sum/redundancySst2Sum) * sigmaSst2(arcNo);  // compute new sigma (for this arc)
      sigmaAccNew(arcNo)  = std::sqrt(ePeAccSum /redundancyAccSum)  * sigmaAcc(arcNo);   // compute new sigma (for this arc)
    }

    // variance component estimation (pod1)
    // ------------------------------------
    if(countPod1)
    {
      std::vector<UInt> index(countPod1);
      for(UInt i=0; i<index.size(); i++)
        index.at(i) = static_cast<UInt>(std::round((observationArc.at(arcNo).timesPod1.at(i)-observationArc.at(arcNo).timesPod1.at(0)).seconds()/samplingPod1));
      Double ePeSum=0, redundancySum=0;
      Matrix R;
      Vector WWe;
      Vce::redundancy(WPod1, We.row(idxPod1, 3*countPod1), WAz.row(idxPod1, 3*countPod1), WB.row(idxPod1, 3*countPod1), R, WWe);
      Vce::psd(R, WWe, index, sigmaPod1(arcNo), CosTransformPod1, PsdPod1, ePePod1, redundancyPod1, ePeSum, redundancySum);
      sigmaPod1New(arcNo) = std::sqrt(ePeSum/redundancySum) * sigmaPod1(arcNo);  // compute new sigma (for this arc)
    }

    // variance component estimation (pod2)
    // ------------------------------------
    if(countPod2)
    {
      std::vector<UInt> index(countPod2);
      for(UInt i=0; i<index.size(); i++)
        index.at(i) = static_cast<UInt>(std::round((observationArc.at(arcNo).timesPod2.at(i)-observationArc.at(arcNo).timesPod2.at(0)).seconds()/samplingPod2));
      Double ePeSum=0, redundancySum=0;
      Matrix R;
      Vector WWe;
      Vce::redundancy(WPod2, We.row(idxPod2, 3*countPod2), WAz.row(idxPod2, 3*countPod2), WB.row(idxPod2, 3*countPod2), R, WWe);
      Vce::psd(R, WWe, index, sigmaPod2(arcNo), CosTransformPod2, PsdPod2, ePePod2, redundancyPod2, ePeSum, redundancySum);
      sigmaPod2New(arcNo) = std::sqrt(ePeSum/redundancySum) * sigmaPod2(arcNo);  // compute new sigma (for this arc)
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void PreprocessingDualSst::varianceComponentEstimationPsd(
                                    const_MatrixSliceRef R, const_MatrixSliceRef WWe, const_MatrixSliceRef CosTransform,
                                    Double sigma1, Double sigma2, Double sigma12,
                                    const_MatrixSliceRef Psd1, const_MatrixSliceRef Psd2, const_MatrixSliceRef Psd12,
                                    MatrixSliceRef ePe1,  MatrixSliceRef redundancy1,  Double &ePe1Sum,  Double &redundancy1Sum,
                                    MatrixSliceRef ePe2,  MatrixSliceRef redundancy2,  Double &ePe2Sum,  Double &redundancy2Sum,
                                    MatrixSliceRef ePe12, MatrixSliceRef redundancy12, Double &ePe12Sum, Double &redundancy12Sum)
{
  try
  {
    const UInt count = WWe.rows()/2;

    Vector e1(CosTransform.rows()), e2(CosTransform.rows()), e12(CosTransform.rows());
    Vector r1(CosTransform.rows()), r2(CosTransform.rows()), r12(CosTransform.rows());
    for(UInt i=0; i<count; i++)
      for(UInt k=0; k<count; k++)
      {
        e1 ((k<i)?(i-k):(k-i)) += WWe(i,0)       * WWe(k,0);
        e2 ((k<i)?(i-k):(k-i)) += WWe(i+count,0) * WWe(k+count,0);
        e12((k<i)?(i-k):(k-i)) += WWe(i,0)       * WWe(k,0)
                               +  WWe(i+count,0) * WWe(k+count,0)
                               +  WWe(i,0)       * WWe(k+count,0)
                               +  WWe(i+count,0) * WWe(k,0);

        r1 ((k<i)?(i-k):(k-i)) += R(std::min(i,k),             std::max(i,k));
        r2 ((k<i)?(i-k):(k-i)) += R(std::min(i+count,k+count), std::max(i+count,k+count));
        r12((k<i)?(i-k):(k-i)) += R(std::min(i,k),             std::max(i,k))
                               +  R(std::min(i+count,k+count), std::max(i+count,k+count))
                               +  R(std::min(i,k+count),       std::max(i,k+count))
                               +  R(std::min(i+count,k),       std::max(i+count,k));
      }

    for(UInt idFreq=0; idFreq<Psd1.rows(); idFreq++)
    {
      const_MatrixSliceRef cov(CosTransform.column(idFreq));
      const Double ePe1Tmp         = std::pow(sigma1,  2) * Psd1 (idFreq, 0) * inner(e1,  cov);
      const Double ePe2Tmp         = std::pow(sigma2,  2) * Psd2 (idFreq, 0) * inner(e2,  cov);
      const Double ePe12Tmp        = std::pow(sigma12, 2) * Psd12(idFreq, 0) * inner(e12, cov);
      const Double redundancy1Tmp  = std::pow(sigma1,  2) * Psd1 (idFreq, 0) * inner(r1,  cov);
      const Double redundancy2Tmp  = std::pow(sigma2,  2) * Psd2 (idFreq, 0) * inner(r2,  cov);
      const Double redundancy12Tmp = std::pow(sigma12, 2) * Psd12(idFreq, 0) * inner(r12, cov);

      ePe1 (idFreq, 0)        += ePe1Tmp;
      ePe2 (idFreq, 0)        += ePe2Tmp;
      ePe12(idFreq, 0)        += ePe12Tmp;
      redundancy1 (idFreq, 0) += redundancy1Tmp;
      redundancy2 (idFreq, 0) += redundancy2Tmp;
      redundancy12(idFreq, 0) += redundancy12Tmp;
      ePe1Sum                 += ePe1Tmp;
      ePe2Sum                 += ePe2Tmp;
      ePe12Sum                += ePe12Tmp;
      redundancy1Sum          += redundancy1Tmp;
      redundancy2Sum          += redundancy2Tmp;
      redundancy12Sum         += redundancy12Tmp;
    } // for(idFreq)
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
    const UInt countSst  = observationArc.at(arcNo).timesSst.size();
    const UInt countPod1 = observationArc.at(arcNo).timesPod1.size();
    const UInt countPod2 = observationArc.at(arcNo).timesPod2.size();

    UInt obsCount = 0;
    const UInt idxSst1  = obsCount; obsCount += countSst;
    const UInt idxSst2  = obsCount; obsCount += countSst;
    const UInt idxPod1  = obsCount; obsCount += 3*countPod1;
    const UInt idxPod2  = obsCount; obsCount += 3*countPod2;

    // Residuals
    // ---------
    Matrix e = observationArc.at(arcNo).l;
    normals.designMatMult(findInterval(arcNo), -1., observationArc.at(arcNo).A, x, e);

    // eliminate arc dependent parameters
    // ----------------------------------
    Matrix CovSst1, CovSst2, CovAcc, WSst, WPod1, WPod2;
    if(observationArc.at(arcNo).B.size())
    {
      Matrix We = e;
      Matrix WB = observationArc.at(arcNo).B;
      decorrelate(arcNo, observationArc.at(arcNo).timesSst.size(),
                  observationArc.at(arcNo).timesPod1.size(), observationArc.at(arcNo).timesPod2.size(),
                  CovSst1, CovSst2, CovAcc, WSst, WPod1, WPod2, {We, WB});

      Vector tau = QR_decomposition(WB);
      QTransMult(WB, tau, We); // transform observations: l:= Q'l
      Matrix y = We.row(0, tau.rows());
      triangularSolve(1., WB.row(0, tau.rows()), y);
      matMult(-1, observationArc.at(arcNo).B, y, e);
    }
    else
      decorrelate(arcNo, observationArc.at(arcNo).timesSst.size(),
                  observationArc.at(arcNo).timesPod1.size(), observationArc.at(arcNo).timesPod2.size(),
                  CovSst1, CovSst2, CovAcc, WSst, WPod1, WPod2, {});

    // collocation for residual separation
    // -----------------------------------
    Matrix WWe = e.row(idxSst1, 2*countSst);
    triangularSolve(1., WSst.trans(), WWe);
    triangularSolve(1., WSst,         WWe);
    const Matrix eSst1 = CovSst1 *  WWe.row(idxSst1, countSst);
    const Matrix eSst2 = CovSst2 *  WWe.row(idxSst2, countSst);
    const Matrix eAcc  = CovAcc  * (WWe.row(idxSst1, countSst) + WWe.row(idxSst2, countSst));

    // create Sst1 arc
    // ---------------
    for(UInt i=0; i<countSst; i++)
    {
      SatelliteTrackingEpoch epoch;
      epoch.time  = observationArc.at(arcNo).timesSst.at(i);
      epoch.range = epoch.rangeRate = epoch.rangeAcceleration = eSst1(i, 0);
      arcListResidualsSst1.at(arcNo).push_back(epoch);
    }

    // create Sst2 arc
    // ---------------
    for(UInt i=0; i<countSst; i++)
    {
      SatelliteTrackingEpoch epoch;
      epoch.time  = observationArc.at(arcNo).timesSst.at(i);
      epoch.range = epoch.rangeRate = epoch.rangeAcceleration = eSst2(i, 0);
      arcListResidualsSst2.at(arcNo).push_back(epoch);
    }

    // create Acc arc
    // --------------
    for(UInt i=0; i<countSst; i++)
    {
      SatelliteTrackingEpoch epoch;
      epoch.time  = observationArc.at(arcNo).timesSst.at(i);
      epoch.range = epoch.rangeRate = epoch.rangeAcceleration = eAcc(i, 0);
      arcListResidualsAcc.at(arcNo).push_back(epoch);
    }

    // create Pod1 arc
    // ---------------
    for(UInt i=0; i<countPod1; i++)
    {
      OrbitEpoch epoch;
      epoch.time     = observationArc.at(arcNo).timesPod1.at(i);
      epoch.position = Vector3d(e(idxPod1+3*i+0,0), e(idxPod1+3*i+1,0), e(idxPod1+3*i+2,0));
      arcListResidualsPod1.at(arcNo).push_back(epoch);
    }

    // create Pod2 arc
    // ---------------
    for(UInt i=0; i<countPod2; i++)
    {
      OrbitEpoch epoch;
      epoch.time     = observationArc.at(arcNo).timesPod2.at(i);
      epoch.position = Vector3d(e(idxPod2+3*i+0,0), e(idxPod2+3*i+1,0), e(idxPod2+3*i+2,0));
      arcListResidualsPod2.at(arcNo).push_back(epoch);
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

    const UInt   countSst      = observationArc.at(arcNo).timesSst.size();
    const Double thresholdSst1 = huber*huber*covFuncSst1(0,1);
    const Double thresholdSst2 = huber*huber*covFuncSst2(0,1);
    const Double thresholdAcc  = huber*huber*covFuncAcc (0,1);
    for(UInt i=0; i<countSst; i++)
      arcListEpochSigmaSst1.at(arcNo).at(i).sigma = arcListEpochSigmaSst2.at(arcNo).at(i).sigma = arcListEpochSigmaAcc.at(arcNo).at(i).sigma = 0.0;

    for(UInt i=0; i<countSst; i++)
    {
      const Double e2Sst1 = std::pow(arcListResidualsSst1.at(arcNo).at(i).rangeRate, 2);
      const Double e2Sst2 = std::pow(arcListResidualsSst2.at(arcNo).at(i).rangeRate, 2);
      const Double e2Acc  = std::pow(arcListResidualsAcc.at(arcNo).at(i).rangeRate,  2);

      if(e2Sst1 > thresholdSst1) arcListEpochSigmaSst1.at(arcNo).at(i).sigma = std::sqrt(e2Sst1-thresholdSst1);
      if(e2Sst2 > thresholdSst2) arcListEpochSigmaSst2.at(arcNo).at(i).sigma = std::sqrt(e2Sst2-thresholdSst2);
      if(e2Acc  > thresholdAcc)  arcListEpochSigmaAcc.at(arcNo).at(i).sigma  = std::sqrt(e2Acc-thresholdAcc);
    }

    const UInt   countPod1     = observationArc.at(arcNo).timesPod1.size();
    const Double thresholdPod1 = huber*huber*(covFuncPod1(0,1)+covFuncPod1(0,2)+covFuncPod1(0,3));
    for(UInt i=0; i<countPod1; i++)
    {
      const Double e2 = std::pow(arcListResidualsPod1.at(arcNo).at(i).position.x(),2)+
                        std::pow(arcListResidualsPod1.at(arcNo).at(i).position.y(),2)+
                        std::pow(arcListResidualsPod1.at(arcNo).at(i).position.z(),2);
      arcListEpochSigmaPod1.at(arcNo).at(i).sigma = 0.;
      if(e2 > thresholdPod1)
        arcListEpochSigmaPod1.at(arcNo).at(i).sigma = std::sqrt((e2-thresholdPod1)/3);
    }

    const UInt   countPod2     = observationArc.at(arcNo).timesPod2.size();
    const Double thresholdPod2 = huber*huber*(covFuncPod2(0,1)+covFuncPod2(0,2)+covFuncPod2(0,3));
    for(UInt i=0; i<countPod2; i++)
    {
      const Double e2 = std::pow(arcListResidualsPod2.at(arcNo).at(i).position.x(),2)+
                        std::pow(arcListResidualsPod2.at(arcNo).at(i).position.y(),2)+
                        std::pow(arcListResidualsPod2.at(arcNo).at(i).position.z(),2);
      arcListEpochSigmaPod2.at(arcNo).at(i).sigma = 0.;
      if(e2 > thresholdPod2)
        arcListEpochSigmaPod2.at(arcNo).at(i).sigma = std::sqrt((e2-thresholdPod2)/3);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
