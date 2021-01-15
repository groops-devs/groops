/***********************************************/
/**
* @file preprocessingSst.cpp
*
* @brief Estimate covariance function / arc weights.
*
* @author Torsten Mayer-Guerr
* @author Matthias Ellmer
* @author Andreas Kvas
* @date 2011-09-23
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This programs processes satellite-to-satellite-tracking (SST) and orbit observations in a GRACE like configuration.
Three different observation groups are considered separatly: SST and POD1/POD2 for the two satellites.
This program works similar to \program{PreprocessingPod}, see there for details. Here only deviations
in the settings are explained.

Precise orbit data (POD) contain often systematic errors additional to stochastic noise. In this case the
variance component estimation fails and assign too much weight to the data. Therefore an additional
\config{downweightPod} factor can be applied to the standard deviation of POD for the next least squares adjustment
in the iteration. This factor should also applied as \config{sigma} in \configClass{observation}{observationType}
for computation of the final solution e.g. with \program{NormalsSolverVCE}.

Daily variations (or other intervals) of the gravity field can be co-estimated together with the static/monthly
mean gravity field. The \configClass{covarianceSignal:parametrizationGravity}{parametrizationGravityType}
describes the daily gravity field parameters and must agree with the parametrization of the covariance matrices.
Furthermore it must be part of
the \configClass{parametrizationGravity}{parametrizationGravityType} selected in \config{observation}.
The intervals are defined in \config{arcList}, created by \program{InstrumentSynchronize}.
These daily estimates must be constrained by providing
\config{inputfileCovariance} matrices in \config{covarianceSignal}. The first file provides the spatial signal
covariances and the second, third and so on files the spatial and temporal covariances between one day and the
next days respectively. The covariance matrices are the same as used in \program{KalmanFilter} and can be
generated with \program{Gravityfield2EmpiricalCovariance}. It is possible to scale the covariance matrices with
\config{covarianceSignal:sigma} and to estimate the factor by means of variance component estimation.
The daily (or in general interval) estimates can be written with \configFile{outputfileSolutionSignal}{matrix}.
For every interval a file is written with an index appended (counting from zero) to the file name.
)";

/***********************************************/

#include "programs/program.h"
#include "parallel/matrixDistributed.h"
#include "files/fileArcList.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"
#include "files/fileParameterName.h"
#include "classes/parametrizationGravity/parametrizationGravity.h"
#include "misc/observation/observationMiscSst.h"
#include "misc/varianceComponentEstimation.h"
#include "misc/kalmanProcessing.h"

/***** CLASS ***********************************/

/** @brief Estimate covariance function / arc weights.
* @ingroup programsGroup */
class PreprocessingSst
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);

private:
  ObservationMiscSstPtr observationMisc;
  std::vector<ObservationMiscSst::Arc> observationArc;
  InstrumentFile fileCovEpochPod1, fileCovEpochPod2;

  std::vector<UInt> arcsInterval;
  std::vector<Time> timesInterval;
  std::vector< std::vector<UInt> > arcIndexN;
  std::vector< std::vector<UInt> > arcIndexA;

  Bool estimateCovarianceFunctionVCE;
  Bool estimateVarianceFactorsSstArcWiseCovariances;
  Bool estimateArcSigmas;
  Bool estimateEpochSigmas;
  Bool estimateResiduals;
  Bool estimateArModelSigma;

  // normal equations
  // ----------------
  std::vector<ParameterName> blockParaName; // name of first parameter in each block
  std::vector<UInt> blockInterval;
  std::vector<UInt> blockIndexHighFrequency;
  std::vector<UInt> blockIndexInterval;
  UInt              blockIndexStatic;
  MatrixDistributed normals;
  Matrix n;        // right hand sides
  Matrix x;        // solution
  Matrix Wz;       // monte carlo vector for redundancy computation
  Double lPl;      // =l'Pl, weighted norm of the observations
  UInt   obsCount; // number of observations

  // covariance
  // ----------
  Vector  sigmaSst, sigmaPod1, sigmaPod2;
  Vector  sigmaSstNew, sigmaPod1New, sigmaPod2New;
  std::vector<ObservationSigmaArc> arcListEpochSigmaSst, arcListEpochSigmaPod1, arcListEpochSigmaPod2;
  std::vector<std::vector<Matrix>> covMatrixSst; // Several independant matrices per arc
  Vector  covMatrixSstSigmas, ePeCovMatrixSst, redundancyCovMatrixSst;
  Matrix  covFuncSst, covFuncPod1, covFuncPod2;
  Matrix  PsdSst, PsdPod1, PsdPod2;
  Matrix  ePeSst, ePePod1, ePePod2;
  Matrix  redundancySst, redundancyPod1, redundancyPod2; // one row for each frequency, one column for each component
  Matrix  CosTransformSst, CosTransformPod1, CosTransformPod2;
  Double  samplingSst, samplingPod1, samplingPod2;

  // residuals
  std::vector<SatelliteTrackingArc> arcListResidualsSst;
  std::vector<OrbitArc> arcListResidualsPod1, arcListResidualsPod2;

  void initNormals(UInt countSignalParameter, const std::vector<ParameterName> &paraNameSignal, Parallel::CommunicatorPtr comm);
  void computeObservationEquation(UInt arcNo);
  void buildNormals(UInt arcNo);
  void computeRedundancies(UInt arcNo);
  void computeResiduals(UInt arcNo);
  void computeEpochSigmas(UInt arcNo);
  Arc  collectArcSst (UInt arcNo) const {return arcListResidualsSst.at(arcNo);}
  Arc  collectArcPod1(UInt arcNo) const {return arcListResidualsPod1.at(arcNo);}
  Arc  collectArcPod2(UInt arcNo) const {return arcListResidualsPod2.at(arcNo);}
  Arc  collectEpochSigmasSst (UInt arcNo) const {return arcListEpochSigmaSst.at(arcNo);}
  Arc  collectEpochSigmasPod1(UInt arcNo) const {return arcListEpochSigmaPod1.at(arcNo);}
  Arc  collectEpochSigmasPod2(UInt arcNo) const {return arcListEpochSigmaPod2.at(arcNo);}

  // Arc-wise SST covariance matrices
  void initArcWiseSSTCovarianceMatrices(UInt arcNo, const std::vector<FileName> names);
  Matrix arcWiseSSTCovarianceMatrix(UInt arcNo) const;
};

GROOPS_REGISTER_PROGRAM(PreprocessingSst, PARALLEL, "Estimate covariance function / arc weights.", Preprocessing, Covariance, Residuals, KalmanFilter)

/***********************************************/
/***********************************************/

void PreprocessingSst::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    estimateCovarianceFunctionVCE = estimateArcSigmas = estimateEpochSigmas = estimateResiduals = FALSE;
    estimateArModelSigma = FALSE;
    estimateVarianceFactorsSstArcWiseCovariances = FALSE;

    FileName fileNameSolution,         fileNameSigmax,            fileNameParaName;
    FileName fileNameHighFrequency,    fileNameInterval;
    FileName fileNameOutArcSigmaSst,   fileNameOutArcSigmaPod1,   fileNameOutArcSigmaPod2;
    FileName fileNameOutEpochSigmaSst, fileNameOutEpochSigmaPod1, fileNameOutEpochSigmaPod2;
    FileName fileNameOutCovSst,        fileNameOutCovPod1,        fileNameOutCovPod2;
    FileName fileNameOutVarianceFactorsSst;
    FileName fileNameResidualsSst,     fileNameResidualsPod1,     fileNameResidualsPod2;

    FileName fileNameArcList;
    FileName fileNameInArcSigmaSst,   fileNameInArcSigmaPod1,   fileNameInArcSigmaPod2;
    FileName fileNameInEpochSigmaSst, fileNameInEpochSigmaPod1, fileNameInEpochSigmaPod2;
    FileName fileNameInCovSst,        fileNameInCovPod1,        fileNameInCovPod2;
    std::vector<FileName> fileNamesCovarianceMatricesSst;
    FileName fileNameInCovMatrixSigmas;

    FileName fileNameInCovEpochPod1,  fileNameInCovEpochPod2;
    Double   sigma0Sst=1, sigma0Pod1=1, sigma0Pod2=1;
    Double   adjustmentThreshold;
    UInt     iterCount;
    Double   downweightPod;

    Double sigma0HighFrequency = 1.0;
    AutoregressiveModelSequencePtr arSequence;
    std::vector<ParameterName> parameterNamesHighFrequency;

    // loops
    std::string iterVariableName;

    renameDeprecatedConfig(config, "arcList", "inputfileArcList", date2time(2020, 7, 7));

    readConfig(config, "outputfileSolution",              fileNameSolution,      Config::OPTIONAL, "", "estimated parameter vector (static part only)");
    readConfig(config, "outputfileSigmax",                fileNameSigmax,        Config::OPTIONAL, "", "standard deviations of the parameters (sqrt of the diagonal of the inverse normal equation)");
    readConfig(config, "outputfileParameterName",         fileNameParaName,      Config::OPTIONAL, "", "estimated signal parameters (index is appended)");
    readConfig(config, "outputfileSolutionHighFrequency", fileNameHighFrequency, Config::OPTIONAL, "", "estimated high-frequency gravity field parameters (index is appended for each interval)");
    readConfig(config, "outputfileSolutionInterval",      fileNameInterval,      Config::OPTIONAL, "", "other interval parameters (index is appended)");
    if(readConfigSequence(config, "estimateArcSigmas", Config::OPTIONAL, "", ""))
    {
      estimateArcSigmas = TRUE;
      readConfig(config, "outputfileSigmasPerArcSst",  fileNameOutArcSigmaSst,  Config::OPTIONAL, "", "accuracies of each arc (SST)");
      readConfig(config, "outputfileSigmasPerArcPod1", fileNameOutArcSigmaPod1, Config::OPTIONAL, "", "accuracies of each arc (POD1)");
      readConfig(config, "outputfileSigmasPerArcPod2", fileNameOutArcSigmaPod2, Config::OPTIONAL, "", "accuracies of each arc (POD2)");
      endSequence(config);
    }
    if(readConfigSequence(config, "estimateEpochSigmas", Config::OPTIONAL, "", ""))
    {
      estimateEpochSigmas = TRUE;
      readConfig(config, "outputfileSigmasPerEpochSst",  fileNameOutEpochSigmaSst,  Config::OPTIONAL, "", "accuracies of each epoch (SST)");
      readConfig(config, "outputfileSigmasPerEpochPod1", fileNameOutEpochSigmaPod1, Config::OPTIONAL, "", "accuracies of each epoch (POD1)");
      readConfig(config, "outputfileSigmasPerEpochPod2", fileNameOutEpochSigmaPod2, Config::OPTIONAL, "", "accuracies of each epoch (POD2)");
      endSequence(config);
    }
    if(readConfigSequence(config, "estimateCovarianceFunctions", Config::OPTIONAL, "", ""))
    {
      estimateCovarianceFunctionVCE = TRUE;
      readConfig(config, "outputfileCovarianceFunctionSst",  fileNameOutCovSst,  Config::OPTIONAL, "", "covariance function");
      readConfig(config, "outputfileCovarianceFunctionPod1", fileNameOutCovPod1, Config::OPTIONAL, "", "covariance functions for along, cross, radial direction");
      readConfig(config, "outputfileCovarianceFunctionPod2", fileNameOutCovPod2, Config::OPTIONAL, "", "covariance functions for along, cross, radial direction");
      endSequence(config);
    }
    if(readConfigSequence(config, "estimateSstArcWiseCovariancesVarianceFactors", Config::OPTIONAL, "", ""))
    {
      estimateVarianceFactorsSstArcWiseCovariances = TRUE;
      readConfig(config, "outputfileVarianceFactors",  fileNameOutVarianceFactorsSst,  Config::OPTIONAL, "", "one variance factor per matrix");
      endSequence(config);
    }
    if(readConfigSequence(config, "computeResiduals", Config::OPTIONAL, "", ""))
    {
      estimateResiduals = TRUE;
      readConfig(config, "outputfileSstResiduals",  fileNameResidualsSst,  Config::OPTIONAL, "", "");
      readConfig(config, "outputfilePod1Residuals", fileNameResidualsPod1, Config::OPTIONAL, "", "");
      readConfig(config, "outputfilePod2Residuals", fileNameResidualsPod2, Config::OPTIONAL, "", "");
      endSequence(config);
    }
    readConfig(config, "observation", observationMisc, Config::MUSTSET,  "", "");
    if(readConfigSequence(config, "covarianceSst", Config::MUSTSET, "", ""))
    {
      readConfig(config, "sigma",                              sigma0Sst,                      Config::DEFAULT,  "1", "apriori factor of covariance function");
      readConfig(config, "inputfileSigmasPerArc",              fileNameInArcSigmaSst,          Config::OPTIONAL, "",  "apriori different accuaries for each arc (multiplicated with sigma)");
      readConfig(config, "inputfileSigmasPerEpoch",            fileNameInEpochSigmaSst,        Config::OPTIONAL, "",  "apriori different accuaries for each epoch");
      readConfig(config, "inputfileCovarianceFunction",        fileNameInCovSst,               Config::OPTIONAL, "",  "approximate covariances in time");
      readConfig(config, "inputfileCovarianceMatrixArc",       fileNamesCovarianceMatricesSst, Config::OPTIONAL, "",  "Must be given per sst arc with correct dimensions.");
      readConfig(config, "inputfileSigmasCovarianceMatrixArc", fileNameInCovMatrixSigmas,      Config::OPTIONAL, "",  "Vector with one sigma for each <inputfileCovarianceMatrixArc>");
      readConfig(config, "sampling",                           samplingSst,                    Config::DEFAULT,  "5", "[seconds] sampling of the covariance function");
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
    if(readConfigSequence(config, "estimateHighFrequencyVariations", Config::OPTIONAL, "", "co-estimate high-frequency gravity field variations"))
    {
      ParametrizationGravityPtr parametrization;

      readConfig(config, "estimateSigma",               estimateArModelSigma, Config::DEFAULT,  "0",   "estimate factor of covariance  function (VCE)");
      readConfig(config, "sigma",                       sigma0HighFrequency,  Config::DEFAULT,  "1.0", "a-priori standard deviation for AR model regularization");
      readConfig(config, "autoregressiveModelSequence", arSequence,           Config::MUSTSET,  "",    "AR model sequence for constraining high frequency gravity variations");
      readConfig(config, "parametrization",              parametrization,       Config::MUSTSET,  "",    "must be set in observation too");
      endSequence(config);

      if(!isCreateSchema(config))
        parametrization->parameterName(parameterNamesHighFrequency);
    }
    readConfig(config, "inputfileArcList",       fileNameArcList,    Config::OPTIONAL, "",  "list to correspond points of time to arc numbers");
    readConfig(config, "adjustmentThreshold",    adjustmentThreshold,Config::DEFAULT,  "0", "Adjustment factor threshold: Iteration will be stopped once both SST and POD adjustment factors are under this threshold");
    readConfig(config, "iterationCount",         iterCount,          Config::DEFAULT,  "3", "(maximum) number of iterations for the estimation of calibration parameter and error PSD");
    readConfig(config, "downweightPod",          downweightPod,      Config::DEFAULT,  "1", "downweight factor for POD");
    readConfig(config, "variableNameIterations", iterVariableName,   Config::OPTIONAL, "",  "All output fileNames in preprocessing iteration are expanded with this variable prior to writing to disk");

    if(isCreateSchema(config)) return;

    // =============================================

    // init
    // ----
    const UInt arcCount        = observationMisc->arcCount();
    const UInt countAParameter = observationMisc->parameterCount();

    arcsInterval = {0, arcCount};
    if(!fileNameArcList.empty())
    {
      logStatus<<"read arc list <"<<fileNameArcList<<">"<<Log::endl;
      readFileArcList(fileNameArcList, arcsInterval, timesInterval);
    }
    else if(arSequence != nullptr)
      throw(Exception("arcList must be given, if autoregressiveModelSequence is set"));

    // init arc sigmas
    // ---------------
    sigmaSst = sigmaPod1 = sigmaPod2 = Vector(arcCount);
    for(UInt arcNo=0; arcNo<arcCount; arcNo++)
      sigmaSst(arcNo) = sigmaPod1(arcNo) = sigmaPod2(arcNo) = 1.0;
    if(!fileNameInArcSigmaSst.empty())  readFileMatrix(fileNameInArcSigmaSst,  sigmaSst);
    if(!fileNameInArcSigmaPod1.empty()) readFileMatrix(fileNameInArcSigmaPod1, sigmaPod1);
    if(!fileNameInArcSigmaPod1.empty()) readFileMatrix(fileNameInArcSigmaPod2, sigmaPod2);
    if(sigmaSst.rows()  != arcCount) throw(Exception("sigmasPerArc (SST) contains wrong number of arcs"));
    if(sigmaPod1.rows() != arcCount) throw(Exception("sigmasPerArc (POD1) contains wrong number of arcs"));
    if(sigmaPod2.rows() != arcCount) throw(Exception("sigmasPerArc (POD2) contains wrong number of arcs"));

    fileCovEpochPod1.open(fileNameInCovEpochPod1);
    fileCovEpochPod2.open(fileNameInCovEpochPod2);

    // ===================================================

    // normal equations of high-frequency variations
    // ---------------------------------------------
    UInt countHighFrequencyParameters = 0;
    std::vector< std::vector< std::vector<Matrix> > > normalsHighFrequency;
    if(arSequence != nullptr)
    {
      logStatus<<"initialize normal equations for high-frequency gravity field parameters"<<Log::endl;
      countHighFrequencyParameters = arSequence->dimension();
      normalsHighFrequency = arSequence->normalEquationSequence();
    }

    // =============================================

    // Init normal equations
    // ---------------------
    logStatus<<"Init normal equations"<<Log::endl;
    initNormals(countHighFrequencyParameters, parameterNamesHighFrequency, comm);

    // ===================================================

    // parameter names
    // ---------------
    if(!fileNameParaName.empty() && Parallel::isMaster(comm))
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
    std::vector<UInt> processNo = Parallel::forEachInterval(arcCount, arcsInterval, [this](UInt arcNo) {computeObservationEquation(arcNo);}, comm);
    observationMisc = ObservationMiscSstPtr(nullptr);

    // =============================================

    // read SST covariance matrices
    // ----------------------------
    const UInt covarianceSstArcWiseCount = fileNamesCovarianceMatricesSst.size();
    covMatrixSst.resize(arcCount);
    covMatrixSstSigmas = Vector(covarianceSstArcWiseCount);
    for(UInt i=0; i<covarianceSstArcWiseCount; i++)
      covMatrixSstSigmas(i) = 1.0;

    if(!fileNameInCovMatrixSigmas.empty())
    {
      readFileMatrix(fileNameInCovMatrixSigmas, covMatrixSstSigmas);
      if(covMatrixSstSigmas.rows() != covarianceSstArcWiseCount)
        throw(Exception("Number of sigmas not compatible with number of given arc-wise SST covariance matrices"));
    }

    if(covarianceSstArcWiseCount)
    {
      logStatus<<"read arc-wise sst covariance matrices"<<Log::endl;
      Parallel::forEachProcess(arcCount, [this, &fileNamesCovarianceMatricesSst](UInt arcNo) {initArcWiseSSTCovarianceMatrices(arcNo, fileNamesCovarianceMatricesSst);}, processNo, comm);
    }

    // =============================================

    // count used blocks
    // -----------------
    logStatus<<"setup normal equations"<<Log::endl;
    Matrix BlockUsed(normals.blockCount(), Matrix::SYMMETRIC);
    for(UInt arcNo=0; arcNo<arcCount; arcNo++)
      for(UInt i=0; i<arcIndexN.at(arcNo).size(); i++)
        for(UInt k=0; k<arcIndexN.at(arcNo).size(); k++)
          BlockUsed(arcIndexN.at(arcNo).at(i), arcIndexN.at(arcNo).at(k)) = 1;
    Parallel::reduceSum(BlockUsed, 0, comm);
    Parallel::broadCast(BlockUsed, 0, comm);

    for(UInt i=0; i<normals.blockCount(); i++)
      for(UInt k=i; k<normals.blockCount(); k++)
        if((BlockUsed(i,k)>0) && (!normals.isBlockUsed(i,k)))
          normals.setBlock(i,k);

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
    Parallel::reduceMax(covLengthSst,  0, comm); Parallel::broadCast(covLengthSst,  0, comm);
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
    covFuncSst  = Vce::readCovarianceFunction(fileNameInCovSst,  covLengthSst,  1, samplingSst);
    covFuncPod1 = Vce::readCovarianceFunction(fileNameInCovPod1, covLengthPod1, 3, samplingPod1);
    covFuncPod2 = Vce::readCovarianceFunction(fileNameInCovPod2, covLengthPod2, 3, samplingPod2);
    covFuncSst.column(1,1)  *= pow(sigma0Sst, 2);
    covFuncPod1.column(1,3) *= pow(sigma0Pod1,2);
    covFuncPod2.column(1,3) *= pow(sigma0Pod2,2);
    PsdSst  = CosTransformSst  * covFuncSst.column(1,1);
    PsdPod1 = CosTransformPod1 * covFuncPod1.column(1,3);
    PsdPod2 = CosTransformPod2 * covFuncPod2.column(1,3);

    // =============================================

    // init epoch sigmas
    // -----------------
    arcListEpochSigmaSst.resize(arcCount);
    arcListEpochSigmaPod1.resize(arcCount);
    arcListEpochSigmaPod2.resize(arcCount);

    if(estimateEpochSigmas)
    {
      InstrumentFile fileSst (fileNameInEpochSigmaSst);
      InstrumentFile filePod1(fileNameInEpochSigmaPod1);
      InstrumentFile filePod2(fileNameInEpochSigmaPod2);

      for(UInt arcNo=0; arcNo<arcCount; arcNo++)
        if(Parallel::myRank(comm) == processNo.at(arcNo))
        {
          arcListEpochSigmaSst.at(arcNo)  = fileSst.readArc(arcNo);
          arcListEpochSigmaPod1.at(arcNo) = filePod1.readArc(arcNo);
          arcListEpochSigmaPod2.at(arcNo) = filePod2.readArc(arcNo);

          if(arcListEpochSigmaSst.at(arcNo).size()==0)
            for(UInt i=0; i<observationArc.at(arcNo).timesSst.size(); i++)
            {
              ObservationSigmaEpoch epoch;
              epoch.time = observationArc.at(arcNo).timesSst.at(i);
              arcListEpochSigmaSst.at(arcNo).push_back(epoch);
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

    Double sigma2_signal = pow(sigma0HighFrequency, 2);
    for(UInt iter=0; iter<iterCount; iter++)
    {
      logInfo<<"starting iteration "<<iter<<Log::endl;
      bool thresholdReached = false;

      if(iterVariableName.size())
        variableIteration[iterVariableName]->setValue(iter);

      // solve normal equations
      // ----------------------
      if(countAParameter)
      {
        logStatus<<"accumulate system of normal equations"<<Log::endl;
        normals.setNull();
        n.setNull();
        lPl      = 0;
        obsCount = 0;
        Parallel::forEachProcess(arcCount, [this](UInt arcNo) {buildNormals(arcNo);}, processNo, comm);

        // collect system of normal equations
        // ----------------------------------
        logStatus<<"collect system of normal equations"<<Log::endl;
        Parallel::reduceSum(n, 0, comm);
        Parallel::reduceSum(obsCount, 0, comm);
        Parallel::reduceSum(lPl, 0, comm);
        normals.reduceSum();

        // =============================================

        // add normals of AR model
        // -----------------------
        if(blockIndexHighFrequency.size())
        {
          logStatus<<"normals of process dynamic"<<Log::endl;
          if(Parallel::isMaster(comm))
            obsCount += countHighFrequencyParameters * blockIndexHighFrequency.size();
          for(UInt id=0; id<blockIndexHighFrequency.size(); id++)
          {
            const UInt idx = std::min(id, normalsHighFrequency.size()-1);
            for(UInt i=0; i<normalsHighFrequency.at(idx).size(); i++)
              for(UInt k=i; k<normalsHighFrequency.at(idx).at(i).size(); k++)
              {
                normals.setBlock(blockIndexHighFrequency.at(id+i-idx), blockIndexHighFrequency.at(id+k-idx));
                if(normals.isMyRank(blockIndexHighFrequency.at(id+i-idx), blockIndexHighFrequency.at(id+k-idx)))
                  axpy(1./sigma2_signal, normalsHighFrequency.at(idx).at(i).at(k), normals.N(blockIndexHighFrequency.at(id+i-idx), blockIndexHighFrequency.at(id+k-idx)));
              }
          }
        } // if(blockIndexSignal.size())

        // =============================================

        // Regularize not used parameters
        // ------------------------------
        UInt countRegul = 0;
        for(UInt i=0; i<normals.blockCount(); i++)
          if(normals.isMyRank(i,i))
          {
            Matrix &N = normals.N(i,i);
            for(UInt k=0; k<N.rows(); k++)
              if(N(k,k) == 0)
              {
                N(k,k) = 1.;
                countRegul++;
              }
          }
        Parallel::reduceSum(countRegul, 0, comm);
        if(Parallel::isMaster(comm) && countRegul)
          logWarning<<countRegul<<" parameters are not used"<<Log::endl;
        Parallel::barrier(comm);

        // =============================================

        // cholesky and forward step
        // -------------------------
        logStatus<<"solve system of normal equations"<<Log::endl;
        x = normals.solve(n, TRUE/*timing*/);
        Parallel::broadCast(x, 0, comm);
        if(Parallel::isMaster(comm))
          logInfo<<"  aposteriori sigma = "<<sqrt((lPl-inner(x, n))/(obsCount-normals.parameterCount()))<<Log::endl;

        // N contains now the cholesky decomposition
        Wz = Vce::monteCarlo(normals.parameterCount(), 100); // monte carlo vector for VCE
        normals.triangularSolve(Wz);
        Parallel::broadCast(Wz, 0, comm);

        if(Parallel::isMaster(comm) && !fileNameSolution.empty())
        {
          logStatus<<"write solution to <"<<fileNameSolution(variableIteration)<<">"<<Log::endl;
          writeFileMatrix(fileNameSolution(variableIteration), x.row(normals.blockIndex(blockIndexStatic), normals.parameterCount()-normals.blockIndex(blockIndexStatic)));
        }

        if(!fileNameSigmax.empty())
        {
          logStatus<<"inverte cholesky matrix and write standard deviations to <"<<fileNameSigmax(variableIteration)<<">"<<Log::endl;
          for(UInt i=blockIndexStatic; i<normals.blockCount(); i++)
            for(UInt k=i; k<normals.blockCount(); k++)
              if(normals.rank(i,k) != 0)
              {
                if(normals.isMyRank(i,k))
                  Parallel::send(normals.N(i,k), 0, comm);
                else if(Parallel::isMaster(comm))
                  Parallel::receive(normals.N(i,k), normals.rank(i,k), comm);
              }
          if(Parallel::isMaster(comm))
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
            writeFileMatrix(fileNameSigmax(variableIteration), diagonal);
          }
          for(UInt i=blockIndexStatic; i<normals.blockCount(); i++)
            for(UInt k=i; k<normals.blockCount(); k++)
              if(!normals.isMyRank(i,k))
                normals.N(i,k) = Matrix();
        } // if(!fileNameSigmax.empty())

        if(Parallel::isMaster(comm) && blockIndexHighFrequency.size() && !fileNameHighFrequency.empty())
        {
          logStatus<<"write signal solution to <"<<fileNameHighFrequency(variableIteration)<<">"<<Log::endl;
          for(UInt i=0; i<blockIndexHighFrequency.size(); i++)
            writeFileMatrix(fileNameHighFrequency(variableIteration).appendBaseName(i%".%02i"s), x.row(normals.blockIndex(blockIndexHighFrequency.at(i)), normals.blockSize(blockIndexHighFrequency.at(i))));
        } // if(!fileNameSignal.empty())

        if(Parallel::isMaster(comm) && blockIndexInterval.size() && !fileNameInterval.empty())
        {
          logStatus<<"write interval solution to <"<<fileNameInterval(variableIteration)<<">"<<Log::endl;
          for(UInt i=0; i<blockIndexInterval.size(); i++)
            writeFileMatrix(fileNameInterval(variableIteration).appendBaseName(i%".%02i"s), x.row(normals.blockIndex(blockIndexInterval.at(i)), normals.blockSize(blockIndexInterval.at(i))));
        } // if(!fileNameInterval.empty())
      } // if(countAParameter)
      Parallel::barrier(comm);

      // =============================================

      if(estimateResiduals || estimateEpochSigmas)
      {
        logStatus<<"compute residuals"<<Log::endl;
        arcListResidualsSst.clear();  arcListResidualsSst.resize(arcCount);
        arcListResidualsPod1.clear(); arcListResidualsPod1.resize(arcCount);
        arcListResidualsPod2.clear(); arcListResidualsPod2.resize(arcCount);
        Parallel::forEachProcess(arcCount, [this](UInt arcNo) {computeResiduals(arcNo);}, processNo, comm);
        Parallel::forEachProcess(arcListResidualsSst,  [this](UInt arcNo) {return collectArcSst (arcNo);}, processNo, comm);
        Parallel::forEachProcess(arcListResidualsPod1, [this](UInt arcNo) {return collectArcPod1(arcNo);}, processNo, comm);
        Parallel::forEachProcess(arcListResidualsPod2, [this](UInt arcNo) {return collectArcPod2(arcNo);}, processNo, comm);

        if(Parallel::isMaster(comm) && (!fileNameResidualsSst.empty()))
        {
          logStatus<<"write residual (SST) file <"<<fileNameResidualsSst(variableIteration)<<">"<<Log::endl;
          InstrumentFile::write(fileNameResidualsSst(variableIteration), arcListResidualsSst);
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

      if(estimateArModelSigma && blockIndexHighFrequency.size())
      {
        logStatus<<"compute variance factor of AR model"<<Log::endl;
        Matrix Nx(x.rows(), 1);
        Matrix NWz(Wz.rows(), Wz.columns());
        logTimerStart;
        for(UInt id=0; id<blockIndexHighFrequency.size(); id++)
        {
          logTimerLoop(id, blockIndexHighFrequency.size());
          const UInt idx = std::min(id, normalsHighFrequency.size()-1);
          for(UInt i=0; i<normalsHighFrequency.at(idx).size(); i++)
            for(UInt k=i; k<normalsHighFrequency.at(idx).at(i).size(); k++)
              if(normals.isMyRank(id+i-idx,id+k-idx))
              {
                matMult(1., normalsHighFrequency.at(idx).at(i).at(k),  x.row(normals.blockIndex(blockIndexHighFrequency.at(id+k-idx)), normals.blockSize(blockIndexHighFrequency.at(id+k-idx))),
                                                            Nx.row(normals.blockIndex(blockIndexHighFrequency.at(id+i-idx)), normals.blockSize(blockIndexHighFrequency.at(id+i-idx))));
                matMult(1., normalsHighFrequency.at(idx).at(i).at(k), Wz.row(normals.blockIndex(blockIndexHighFrequency.at(id+k-idx)), normals.blockSize(blockIndexHighFrequency.at(id+k-idx))),
                                                          NWz.row(normals.blockIndex(blockIndexHighFrequency.at(id+i-idx)), normals.blockSize(blockIndexHighFrequency.at(id+i-idx))));
                if(k>i) // extend symmetric
                {
                  matMult(1., normalsHighFrequency.at(idx).at(i).at(k).trans(),  x.row(normals.blockIndex(blockIndexHighFrequency.at(id+i-idx)), normals.blockSize(blockIndexHighFrequency.at(id+i-idx))),
                                                                      Nx.row(normals.blockIndex(blockIndexHighFrequency.at(id+k-idx)), normals.blockSize(blockIndexHighFrequency.at(id+k-idx))));
                  matMult(1., normalsHighFrequency.at(idx).at(i).at(k).trans(), Wz.row(normals.blockIndex(blockIndexHighFrequency.at(id+i-idx)), normals.blockSize(blockIndexHighFrequency.at(id+i-idx))),
                                                                    NWz.row(normals.blockIndex(blockIndexHighFrequency.at(id+k-idx)), normals.blockSize(blockIndexHighFrequency.at(id+k-idx))));
                }
              }
        } // for(id)
        Parallel::reduceSum(Nx,  0, comm);
        Parallel::reduceSum(NWz, 0, comm);
        logTimerLoopEnd(blockIndexHighFrequency.size());

        if(Parallel::isMaster(comm))
        {
          const Double ePe = inner(x,Nx);
          const Double r   = countHighFrequencyParameters*blockIndexHighFrequency.size() - 1./sigma2_signal*inner(Wz,NWz);
          const Double s2  = ePe/r;
          logInfo<<"  Process dynamic sigma: sqrt("<<ePe<<"/("<<countHighFrequencyParameters*blockIndexHighFrequency.size()<<"-"<<1./sigma2_signal*inner(Wz,NWz)<<"))="<<sqrt(s2)<<Log::endl;
          if((s2==s2) && (s2>0))
            sigma2_signal = s2;
        }
        Parallel::broadCast(sigma2_signal, 0, comm);
      } // if(blockIndexSignal.size())

      // =============================================

      // compute redundancies
      // --------------------
      if((estimateArcSigmas || estimateCovarianceFunctionVCE || estimateVarianceFactorsSstArcWiseCovariances))
      {
        logStatus<<"compute redundancies"<<Log::endl;
        sigmaSstNew = sigmaPod1New = sigmaPod2New = Vector(arcCount);
        ePeSst  = redundancySst  = Matrix(covLengthSst,  1);
        ePePod1 = redundancyPod1 = Matrix(covLengthPod1, 3); // for x,y,z
        ePePod2 = redundancyPod2 = Matrix(covLengthPod2, 3); // for x,y,z
        ePeCovMatrixSst = redundancyCovMatrixSst = Vector(covarianceSstArcWiseCount);
        Parallel::forEachProcess(arcCount, [this](UInt arcNo) {computeRedundancies(arcNo);}, processNo, comm);
      }

      // =============================================

      // sigmas per arc
      // --------------
      if(estimateArcSigmas)
      {
        Parallel::reduceSum(sigmaSstNew,  0, comm);
        Parallel::reduceSum(sigmaPod1New, 0, comm);
        Parallel::reduceSum(sigmaPod2New, 0, comm);
        if(Parallel::isMaster(comm))
        {
          sigmaSst  = sigmaSstNew;
          sigmaPod1 = sigmaPod1New;
          sigmaPod2 = sigmaPod2New;

          Double sigma0Sst  = Vce::meanSigma(sigmaSst);
          Double sigma0Pod1 = Vce::meanSigma(sigmaPod1);
          Double sigma0Pod2 = Vce::meanSigma(sigmaPod2);

          sigmaSst  *= 1./sigma0Sst;
          sigmaPod1 *= 1./sigma0Pod1;
          sigmaPod2 *= 1./sigma0Pod2;

          logInfo<<"  SST  sigma per arc (median): "<<sigma0Sst<<Log::endl;
          logInfo<<"  POD1 sigma per arc (median): "<<sigma0Pod1<<Log::endl;
          logInfo<<"  POD2 sigma per arc (median): "<<sigma0Pod2<<Log::endl;

          if(!fileNameOutArcSigmaSst.empty())
          {
            logStatus<<"write arc sigma file <"<<fileNameOutArcSigmaSst(variableIteration)<<">"<<Log::endl;
            writeFileMatrix(fileNameOutArcSigmaSst(variableIteration), sigmaSst);
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
        Parallel::broadCast(sigmaSst,  0, comm);
        Parallel::broadCast(sigmaPod1, 0, comm);
        Parallel::broadCast(sigmaPod2, 0, comm);
      } // if(estimateArcSigmas)

      // =============================================

      if(estimateEpochSigmas)
      {
        logStatus<<"compute epoch sigmas"<<Log::endl;
        Parallel::forEachProcess(arcCount, [this](UInt arcNo) {computeEpochSigmas(arcNo);}, processNo, comm);

        logStatus<<"collect epoch sigmas"<<Log::endl;
        Parallel::forEachProcess(arcListEpochSigmaSst,  [this](UInt arcNo) {return collectEpochSigmasSst (arcNo);}, processNo, comm);
        Parallel::forEachProcess(arcListEpochSigmaPod1, [this](UInt arcNo) {return collectEpochSigmasPod1(arcNo);}, processNo, comm);
        Parallel::forEachProcess(arcListEpochSigmaPod2, [this](UInt arcNo) {return collectEpochSigmasPod2(arcNo);}, processNo, comm);

        if(Parallel::isMaster(comm) && (!fileNameOutEpochSigmaSst.empty()))
        {
          logStatus<<"write epoch sigma (SST) file <"<<fileNameOutEpochSigmaSst(variableIteration)<<">"<<Log::endl;
          InstrumentFile::write(fileNameOutEpochSigmaSst(variableIteration), arcListEpochSigmaSst);
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
        Parallel::reduceSum(ePeSst,         0, comm);
        Parallel::reduceSum(ePePod1,        0, comm);
        Parallel::reduceSum(ePePod2,        0, comm);
        Parallel::reduceSum(redundancySst,  0, comm);
        Parallel::reduceSum(redundancyPod1, 0, comm);
        Parallel::reduceSum(redundancyPod2, 0, comm);

        logStatus<<"compute psd through variance component estimation"<<Log::endl;
        if(Parallel::isMaster(comm))
        {
          Double maxFactorSst = 0;
          Double maxFactorPod = 0;
          Vce::estimatePsd(ePeSst,  redundancySst,  PsdSst,  maxFactorSst);
          Vce::estimatePsd(ePePod1, redundancyPod1, PsdPod1, maxFactorPod);
          Vce::estimatePsd(ePePod2, redundancyPod2, PsdPod2, maxFactorPod);

          maxFactorPod /= downweightPod;
          logInfo<<"  max. PSD adjustment factor (SST): "<<maxFactorSst<<Log::endl;
          logInfo<<"  max. PSD adjustment factor (POD): "<<maxFactorPod<<Log::endl;

          if ((maxFactorSst < adjustmentThreshold) && (maxFactorPod < adjustmentThreshold))
          {
            logStatus<<"  adjustment threshold "<<adjustmentThreshold<<" reached after iteration "<<iter+1<<"."<<Log::endl;
            thresholdReached = true;
          }

        } // if(Parallel::isMaster(comm))
        Parallel::broadCast(PsdSst,  0, comm);
        Parallel::broadCast(PsdPod1, 0, comm);
        Parallel::broadCast(PsdPod2, 0, comm);
        Parallel::broadCast(thresholdReached, 0, comm);
        // compute new covariance function
        copy(CosTransformSst  * PsdSst,  covFuncSst.column(1,1));
        copy(CosTransformPod1 * PsdPod1, covFuncPod1.column(1,3));
        copy(CosTransformPod2 * PsdPod2, covFuncPod2.column(1,3));
      } // if(estimateCovarianceFunctions)

      // =============================================

      // estimate variance factor for arc-wise SST covariance matrices
      // -------------------------------------------------------------
      if(estimateVarianceFactorsSstArcWiseCovariances)
      {
        Parallel::reduceSum(ePeCovMatrixSst, 0, comm);
        Parallel::reduceSum(redundancyCovMatrixSst, 0, comm);
        if(Parallel::isMaster(comm))
        {
          for(UInt i=0; i<covarianceSstArcWiseCount; i++)
          {
            Double alpha = std::sqrt(ePeCovMatrixSst(i) / redundancyCovMatrixSst(i));
            covMatrixSstSigmas(i) *= alpha;
            logStatus<<" SST arc-wise variance factor #"<<i<<" (current/total): "<<alpha<<"/"<<covMatrixSstSigmas(i)<<Log::endl;
          }
        }
        Parallel::broadCast(covMatrixSstSigmas, 0, comm);
      }

      // =============================================

      // Write covariance function to file
      // ---------------------------------
      if(Parallel::isMaster(comm) && !fileNameOutCovSst.empty())
      {
        logStatus<<"write covariance function file <"<<fileNameOutCovSst(variableIteration)<<">"<<Log::endl;
        writeFileMatrix(fileNameOutCovSst(variableIteration), covFuncSst);
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

      // Write variance factor to file
      // -----------------------------
      if(Parallel::isMaster(comm) && !fileNameOutVarianceFactorsSst.empty())
      {
        logStatus<<"write arc-wise SST variance factors <"<<fileNameOutVarianceFactorsSst(variableIteration)<<">"<<Log::endl;
        writeFileMatrix(fileNameOutVarianceFactorsSst(variableIteration), covMatrixSstSigmas);
      }

      // bail if the iteration threshold has been reached
      Parallel::broadCast(thresholdReached, 0, comm);
      if(thresholdReached)
        break;

      if((!estimateArcSigmas) && (!estimateCovarianceFunctionVCE))
        break; // iterations not needed

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

void PreprocessingSst::initArcWiseSSTCovarianceMatrices(UInt arcNo, const std::vector<FileName> names)
{
  try
  {
    for(const auto &fileName : names)
    {
      Matrix c;
      readFileMatrix(fileName.appendBaseName(".arc"+arcNo%"%03i"s), c);
      covMatrixSst.at(arcNo).push_back(c);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix PreprocessingSst::arcWiseSSTCovarianceMatrix(UInt arcNo) const
{
  try
  {
    if(!covMatrixSst.at(arcNo).size())
      return Matrix();

    Matrix W = std::pow(covMatrixSstSigmas.at(0),2) * covMatrixSst.at(arcNo).at(0);
    for(UInt i=1; i<covMatrixSst.at(arcNo).size(); i++)
      axpy(std::pow(covMatrixSstSigmas.at(i),2), covMatrixSst.at(arcNo).at(i), W);

    return W;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

 /***********************************************/

void PreprocessingSst::initNormals(UInt countSignalParameter, const std::vector<ParameterName> &paraNameSignal, Parallel::CommunicatorPtr comm)
{
  try
  {
    // parameter names
    // ---------------
    std::vector<ParameterName> paraNameAll, paraNameStatic;
    observationMisc->parameterName(paraNameAll);
    observationMisc->setInterval(date2time(9999,1,1), date2time(9999,1,2));
    observationMisc->parameterName(paraNameStatic);
    std::vector<Bool> unused(paraNameAll.size(), TRUE);

    // static parameters
    // -----------------
    std::vector<Bool> isStatic(paraNameAll.size(), FALSE);
    UInt k=0;
    for(UInt i=0; i<paraNameStatic.size(); i++)
      for(; k<paraNameAll.size(); k++)
        if(paraNameStatic.at(i) == paraNameAll.at(k))
        {
          isStatic.at(k++) = TRUE;
          break;
        }

    // signal parameters
    // -----------------
    std::vector<UInt> indexSignal(paraNameSignal.size(), NULLINDEX);
    k=0;
    for(UInt i=0; i<paraNameSignal.size(); i++)
      for(; k<paraNameAll.size(); k++)
        if(paraNameSignal.at(i) == paraNameAll.at(k))
        {
          indexSignal.at(i) = k;
          break;
        }
    Bool isSignalPartOfStatic = FALSE;
    for(UInt i=0; i<paraNameSignal.size(); i++)
    {
      if(indexSignal.at(i) == NULLINDEX)
        throw(Exception("Signal covariance parameter '"+paraNameSignal.at(i).str()+"' must be set in observation too."));
      if(isStatic.at(indexSignal.at(i)))
        isSignalPartOfStatic = TRUE;
    }

    // examine intervals
    // -----------------
    std::vector<UInt> blockIndex(1, 0);
    for(UInt idInterval=0; idInterval+1<timesInterval.size(); idInterval++)
    {
      std::vector<ParameterName> paraNameInterval;
      observationMisc->setInterval(timesInterval.at(idInterval), timesInterval.at(idInterval+1));
      observationMisc->parameterName(paraNameInterval);

      // interval parameters
      // -------------------
      std::vector<UInt> indexInterval(paraNameInterval.size(), NULLINDEX);
      std::vector<Bool> inInterval(paraNameAll.size(), FALSE);
      UInt k=0;
      for(UInt i=0; i<paraNameInterval.size(); i++)
        for(; k<paraNameAll.size(); k++)
          if(paraNameInterval.at(i) == paraNameAll.at(k))
          {
            indexInterval.at(i) = k;
            inInterval.at(k) = TRUE;
            break;
          }

      // signal parameters
      // -----------------
      for(UInt i=0; i<indexSignal.size(); i++)
        if(inInterval.at(indexSignal.at(i)) && unused.at(indexSignal.at(i)))
        {
          if(!isSignalPartOfStatic)
            for(UInt count=0; count<countSignalParameter; count++)
              unused.at(indexSignal.at(i+count)) = FALSE;
          blockIndexHighFrequency.push_back(blockIndex.size()-1);
          blockParaName.push_back(paraNameAll.at(indexSignal.at(i)));
          blockIndex.push_back( blockIndex.back()+countSignalParameter );
          blockInterval.push_back((isSignalPartOfStatic) ? idInterval : NULLINDEX);
          i += countSignalParameter-1;
        }

      // other interval parameters
      // -------------------------
      for(UInt i=0; i<indexInterval.size(); i++)
        if(unused.at(indexInterval.at(i)) && !isStatic.at(indexInterval.at(i)))
        {
          UInt count = 0;
          while((i+count<indexInterval.size()) && unused.at(indexInterval.at(i+count)) && !isStatic.at(indexInterval.at(i+count)))
            unused.at(indexInterval.at(i+count++)) = FALSE;
          blockIndexInterval.push_back(blockIndex.size()-1);
          blockParaName.push_back(paraNameAll.at(indexInterval.at(i)));
          blockIndex.push_back( blockIndex.back()+count );
          blockInterval.push_back(NULLINDEX);
          i += count-1;
        }
    } // for(idInterval)

    // static part
    // -----------
    blockIndexStatic = blockIndex.size()-1;
    for(UInt i=0; i<isStatic.size(); i++)
      if(isStatic.at(i))
      {
        UInt count = 0;
        while((i+count<isStatic.size()) && isStatic.at(i+count))
          unused.at(i+count++) = FALSE;
        blockParaName.push_back(paraNameAll.at(i));
        blockIndex.push_back( blockIndex.back()+count );
        blockInterval.push_back(NULLINDEX);
        i += count-1;
      }

    // Init normal equations
    // ---------------------
    normals.initEmpty(blockIndex, comm);
    for(UInt i=0; i<normals.blockCount(); i++)
      normals.setBlock(i,i);
    n = Matrix(normals.parameterCount(), observationMisc->rightSideCount());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void PreprocessingSst::computeObservationEquation(UInt arcNo)
{
  try
  {
    // restrict time interval
    // ----------------------
    UInt idxInterval = 0;
    for(; idxInterval<arcsInterval.size(); idxInterval++)
      if(arcNo<arcsInterval.at(idxInterval+1))
        break;

    if(timesInterval.size())
      observationMisc->setInterval(timesInterval.at(idxInterval), timesInterval.at(idxInterval+1));

    // set block index of normal equations
    // -----------------------------------
    std::vector<ParameterName> paraName;
    observationMisc->parameterName(paraName);
    arcIndexN.at(arcNo).clear();
    arcIndexA.at(arcNo).clear();
    for(UInt idxBlock=0; idxBlock<blockParaName.size(); idxBlock++)
      if((blockInterval.at(idxBlock) == idxInterval) || (blockInterval.at(idxBlock) == NULLINDEX))
        for(UInt idxA=0; idxA<paraName.size(); idxA++)
          if(paraName.at(idxA) == blockParaName.at(idxBlock))
          {
            arcIndexN.at(arcNo).push_back( idxBlock );
            arcIndexA.at(arcNo).push_back( idxA );
            break;
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

void PreprocessingSst::buildNormals(UInt arcNo)
{
  try
  {
    if(observationArc.at(arcNo).A.size() == 0)
      return;

    // count observations and calculate index
    // --------------------------------------
    const UInt countSst  = observationArc.at(arcNo).timesSst.size();
    const UInt countPod1 = observationArc.at(arcNo).timesPod1.size();
    const UInt countPod2 = observationArc.at(arcNo).timesPod2.size();

    UInt obsCount = 0;
    const UInt idxSst  = obsCount; obsCount += countSst;
    const UInt idxPod1 = obsCount; obsCount += 3*countPod1;
    const UInt idxPod2 = obsCount; obsCount += 3*countPod2;

    // Decorrelation
    // -------------
    Matrix Wl   = observationArc.at(arcNo).l;
    Matrix WA   = observationArc.at(arcNo).A;
    Matrix WB   = observationArc.at(arcNo).B;
    Matrix WSst = arcWiseSSTCovarianceMatrix(arcNo);
    Matrix WPod;
    Covariance3dArc covEpochPod1 = fileCovEpochPod1.readArc(arcNo);
    Covariance3dArc covEpochPod2 = fileCovEpochPod2.readArc(arcNo);
    if(countSst)  CovarianceSst::decorrelate(observationArc.at(arcNo).timesSst, sigmaSst(arcNo),  arcListEpochSigmaSst.at(arcNo),  covFuncSst,  WSst, {Wl.row(idxSst, countSst),     WA.row(idxSst, countSst),     WB.row(idxSst, countSst)});
    if(countPod1) CovariancePod::decorrelate(observationArc.at(arcNo).pod1,     sigmaPod1(arcNo), arcListEpochSigmaPod1.at(arcNo), covEpochPod1, covFuncPod1, WPod, {Wl.row(idxPod1, 3*countPod1), WA.row(idxPod1, 3*countPod1), WB.row(idxPod1, 3*countPod1)});
    if(countPod2) CovariancePod::decorrelate(observationArc.at(arcNo).pod2,     sigmaPod2(arcNo), arcListEpochSigmaPod2.at(arcNo), covEpochPod2, covFuncPod2, WPod, {Wl.row(idxPod2, 3*countPod2), WA.row(idxPod2, 3*countPod2), WB.row(idxPod2, 3*countPod2)});

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

void PreprocessingSst::computeRedundancies(UInt arcNo)
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
    const UInt idxSst   = obsCount; obsCount += countSst;
    const UInt idxPod1  = obsCount; obsCount += 3*countPod1;
    const UInt idxPod2  = obsCount; obsCount += 3*countPod2;

    // Decorrelation
    // -------------
    Matrix Wl   = observationArc.at(arcNo).l;
    Matrix WA   = observationArc.at(arcNo).A;
    Matrix WB   = observationArc.at(arcNo).B;
    Matrix WSst = arcWiseSSTCovarianceMatrix(arcNo);
    Matrix WPod1, WPod2;
    Covariance3dArc covEpochPod1 = fileCovEpochPod1.readArc(arcNo);
    Covariance3dArc covEpochPod2 = fileCovEpochPod2.readArc(arcNo);
    if(countSst)  CovarianceSst::decorrelate(observationArc.at(arcNo).timesSst, sigmaSst(arcNo),  arcListEpochSigmaSst.at(arcNo),  covFuncSst,  WSst, {Wl.row(idxSst, countSst),     WA.row(idxSst, countSst),     WB.row(idxSst, countSst)});
    if(countPod1) CovariancePod::decorrelate(observationArc.at(arcNo).pod1,     sigmaPod1(arcNo), arcListEpochSigmaPod1.at(arcNo), covEpochPod1, covFuncPod1, WPod1, {Wl.row(idxPod1, 3*countPod1), WA.row(idxPod1, 3*countPod1), WB.row(idxPod1, 3*countPod1)});
    if(countPod2) CovariancePod::decorrelate(observationArc.at(arcNo).pod2,     sigmaPod2(arcNo), arcListEpochSigmaPod2.at(arcNo), covEpochPod2, covFuncPod2, WPod2, {Wl.row(idxPod2, 3*countPod2), WA.row(idxPod2, 3*countPod2), WB.row(idxPod2, 3*countPod2)});

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
      if(countSst)
      {
        const Double redundancy = countSst - quadsum(WAz.row(idxSst, countSst)) - quadsum(WB.row(idxSst, countSst));
        sigmaSstNew(arcNo) = sqrt(quadsum(We.row(idxSst,  countSst))/redundancy) * sigmaSst(arcNo);
      }
      if(countPod1)
      {
        const Double redundancy = 3*countPod1 - quadsum(WAz.row(idxPod1, 3*countPod1)) - quadsum(WB.row(idxPod1, 3*countPod1));
        sigmaPod1New(arcNo) = sqrt(quadsum(We.row(idxPod1,  3*countPod1))/redundancy) * sigmaPod1(arcNo);
      }
      if(countPod2)
      {
        const Double redundancy = 3*countPod2 - quadsum(WAz.row(idxPod2, 3*countPod2)) - quadsum(WB.row(idxPod2, 3*countPod2));
        sigmaPod2New(arcNo) = sqrt(quadsum(We.row(idxPod2,  3*countPod2))/redundancy) * sigmaPod2(arcNo);
      }
      return;
    } // if(!estimateCovarianceFunctionVCE)

    // ============================================

    // variance component estimation (Sst)
    // -----------------------------------
    if(countSst)
    {
      std::vector<UInt> index(countSst);
      for(UInt i=0; i<index.size(); i++)
        index.at(i) = i;
      Double ePeSum=0, redundancySum=0;
      Matrix R;
      Vector WWe;
      Vce::redundancy(WSst, We.row(idxSst,  countSst), WAz.row(idxSst, countSst), WB.row(idxSst, countSst), R, WWe);

      // Toeplitz covariance function
      Vce::psd(R, WWe, index, sigmaSst(arcNo), CosTransformSst, PsdSst, ePeSst, redundancySst, ePeSum, redundancySum);
      sigmaSstNew(arcNo) = std::sqrt(ePeSum/redundancySum) * sigmaSst(arcNo);  // compute new sigma (for this arc)

      // Arc-wise covariance matrices
      if(estimateVarianceFactorsSstArcWiseCovariances)
        for(UInt i=0; i<covMatrixSst.at(arcNo).size(); i++)
          Vce::matrix(R, WWe, std::pow(covMatrixSstSigmas(i),2) * covMatrixSst.at(arcNo).at(i), ePeCovMatrixSst(i), redundancyCovMatrixSst(i));
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

void PreprocessingSst::computeResiduals(UInt arcNo)
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
    const UInt idxSst   = obsCount; obsCount += countSst;
    const UInt idxPod1  = obsCount; obsCount += 3*countPod1;
    const UInt idxPod2  = obsCount; obsCount += 3*countPod2;

    // Residuals
    // ---------
    Matrix e = observationArc.at(arcNo).l;;
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
      Matrix WSst = arcWiseSSTCovarianceMatrix(arcNo);
      Matrix WPod;

      Covariance3dArc covEpochPod1 = fileCovEpochPod1.readArc(arcNo);
      Covariance3dArc covEpochPod2 = fileCovEpochPod2.readArc(arcNo);
      if(countSst)  CovarianceSst::decorrelate(observationArc.at(arcNo).timesSst, sigmaSst(arcNo),  arcListEpochSigmaSst.at(arcNo),  covFuncSst,  WSst, {We.row(idxSst, countSst), WB.row(idxSst, countSst)});
      if(countPod1) CovariancePod::decorrelate(observationArc.at(arcNo).pod1,     sigmaPod1(arcNo), arcListEpochSigmaPod1.at(arcNo), covEpochPod1, covFuncPod1, WPod, {We.row(idxPod1, 3*countPod1), WB.row(idxPod1, 3*countPod1)});
      if(countPod2) CovariancePod::decorrelate(observationArc.at(arcNo).pod2,     sigmaPod2(arcNo), arcListEpochSigmaPod2.at(arcNo), covEpochPod2, covFuncPod2, WPod, {We.row(idxPod2, 3*countPod2), WB.row(idxPod2, 3*countPod2)});

      Vector tau = QR_decomposition(WB);
      QTransMult(WB, tau, We); // transform observations: l:= Q'l
      Matrix y = We.row(0, tau.rows());
      triangularSolve(1., WB.row(0, tau.rows()), y);
      matMult(-1, observationArc.at(arcNo).B, y, e);
    }

    // create Sst arc
    // --------------
    for(UInt i=0; i<countSst; i++)
    {
      SatelliteTrackingEpoch epoch;
      epoch.time  = observationArc.at(arcNo).timesSst.at(i);
      epoch.range = epoch.rangeRate = epoch.rangeAcceleration = e(idxSst+i,0);
      arcListResidualsSst.at(arcNo).push_back(epoch);
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

void PreprocessingSst::computeEpochSigmas(UInt arcNo)
{
  try
  {
    const Double huber = 2.5;

    const UInt   countSst = observationArc.at(arcNo).timesSst.size();
    const Double thresholdSst = huber*huber*covFuncSst(0,1);
    for(UInt i=0; i<countSst; i++)
    {
      const Double e2 = pow(arcListResidualsSst.at(arcNo).at(i).rangeRate, 2);
      arcListEpochSigmaSst.at(arcNo).at(i).sigma = 0.;
      if(e2 > thresholdSst)
        arcListEpochSigmaSst.at(arcNo).at(i).sigma = sqrt(e2-thresholdSst);
    }

    const UInt   countPod1     = observationArc.at(arcNo).timesPod1.size();
    const Double thresholdPod1 = huber*huber*(covFuncPod1(0,1)+covFuncPod1(0,2)+covFuncPod1(0,3));
    for(UInt i=0; i<countPod1; i++)
    {
      const Double e2 = pow(arcListResidualsPod1.at(arcNo).at(i).position.x(),2)+
                        pow(arcListResidualsPod1.at(arcNo).at(i).position.y(),2)+
                        pow(arcListResidualsPod1.at(arcNo).at(i).position.z(),2);
      arcListEpochSigmaPod1.at(arcNo).at(i).sigma = 0.;
      if(e2 > thresholdPod1)
        arcListEpochSigmaPod1.at(arcNo).at(i).sigma = sqrt((e2-thresholdPod1)/3);
    }

    const UInt   countPod2     = observationArc.at(arcNo).timesPod2.size();
    const Double thresholdPod2 = huber*huber*(covFuncPod2(0,1)+covFuncPod2(0,2)+covFuncPod2(0,3));
    for(UInt i=0; i<countPod2; i++)
    {
      const Double e2 = pow(arcListResidualsPod2.at(arcNo).at(i).position.x(),2)+
                        pow(arcListResidualsPod2.at(arcNo).at(i).position.y(),2)+
                        pow(arcListResidualsPod2.at(arcNo).at(i).position.z(),2);
      arcListEpochSigmaPod2.at(arcNo).at(i).sigma = 0.;
      if(e2 > thresholdPod2)
        arcListEpochSigmaPod2.at(arcNo).at(i).sigma = sqrt((e2-thresholdPod2)/3);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/
