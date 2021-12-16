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
This program processes satellite-to-satellite-tracking (SST) and kinematic orbit observations in a GRACE like configuration.
Three different observation groups are considered separately: SST and POD1/POD2 for the two satellites.
This program works similar to \program{PreprocessingPod}, see there for details. Here only deviations
in the settings are explained.

Precise orbit data (POD) often contains systematic errors in addition to stochastic noise. In this case the
variance component estimation fails and assigns too much weight to the POD data. Therefore an additional
\config{downweightPod} factor can be applied to the standard deviation of POD for the next least squares adjustment
in the iteration. This factor should also applied as \config{sigma} in \configClass{observation}{observationType}
for computation of the final solution e.g. with \program{NormalsSolverVCE}.

Short time variations of the gravity field can be co-estimated together with the static/monthly
mean gravity field. The short time parameters must also be set in \configClass{observation:parametrizationGravity}{parametrizationGravityType} and
can then be selected by \configClass{estimateShortTimeVariations:parameterSelection}{parameterSelectorType}.
If these parameters are not time variable, for example when a range of static parameters is selected,
they are set up as constant for each time interval defined in \config{inputfileArcList}. The parameters are constrained by an
\configClass{estimateShortTimeVariations:autoregressiveModelSequence}{autoregressiveModelSequenceType}. The weight of
the constrain equations in terms of the standard deviation can be estimated by means of
Variance Component Estimation (VCE) if \config{estimateShortTimeVariations:estimateSigma} is set.
The mathematical background of this co-estimation can be found in:

Kvas, A., Mayer-Gürr, T. GRACE gravity field recovery with background model uncertainties.
J Geod 93, 2543–2552 (2019). \url{https://doi.org/10.1007/s00190-019-01314-1}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileArcList.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"
#include "files/fileParameterName.h"
#include "classes/parameterSelector/parameterSelector.h"
#include "misc/observation/observationMiscSst.h"
#include "misc/varianceComponentEstimation.h"
#include "misc/kalmanProcessing.h"
#include "misc/normalsShortTimeStaticLongTime.h"

/***** CLASS ***********************************/

/** @brief Estimate covariance function / arc weights.
* @ingroup programsGroup */
class PreprocessingSst
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);

private:
  ObservationMiscSstPtr                observationMisc;
  std::vector<ObservationMiscSst::Arc> observationArc;
  InstrumentFile                       fileCovEpochPod1, fileCovEpochPod2;
  std::vector<UInt>                    arcsInterval;
  std::vector<Time>                    timesInterval;

  Bool estimateCovarianceFunctionVCE;
  Bool estimateSigmasCovSst;
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
  Vector  sigmaSst, sigmaPod1, sigmaPod2;
  Vector  sigmaSstNew, sigmaPod1New, sigmaPod2New;
  std::vector<ObservationSigmaArc> arcListEpochSigmaSst, arcListEpochSigmaPod1, arcListEpochSigmaPod2;
  std::vector<std::vector<Matrix>> CovSst; // Several independant matrices per arc
  Vector  sigmasCovSst, ePeCovSst, redundancyCovSst;
  Matrix  covFuncSst, covFuncPod1, covFuncPod2;
  Matrix  PsdSst, PsdPod1, PsdPod2;
  Matrix  ePeSst, ePePod1, ePePod2;
  Matrix  redundancySst, redundancyPod1, redundancyPod2; // one row for each frequency, one column for each component
  Matrix  CosTransformSst, CosTransformPod1, CosTransformPod2;
  Double  samplingSst, samplingPod1, samplingPod2;

  // residuals
  std::vector<SatelliteTrackingArc> arcListResidualsSst;
  std::vector<OrbitArc> arcListResidualsPod1, arcListResidualsPod2;

  UInt findInterval(UInt arcNo) const;
  void decorrelate(UInt arcNo, UInt countSst, UInt countPod1, UInt countPod2, Matrix &WSst, Matrix &WPod1, Matrix &WPod2, const std::list<MatrixSlice> &A);
  void buildNormals(UInt arcNo);
  void computeRedundancies(UInt arcNo);
  void computeResiduals(UInt arcNo);
  void computeEpochSigmas(UInt arcNo);
};

GROOPS_REGISTER_PROGRAM(PreprocessingSst, PARALLEL, "Estimate covariance function / arc weights.", Preprocessing, Covariance, Residuals, KalmanFilter)

/***********************************************/
/***********************************************/

void PreprocessingSst::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    estimateCovarianceFunctionVCE = estimateArcSigmas = estimateEpochSigmas = estimateResiduals = FALSE;
    estimateSigmaShortTimeModel = FALSE;
    estimateSigmasCovSst = FALSE;

    FileName fileNameSolution,         fileNameSigmax,            fileNameParaName;
    FileName fileNameOutArcSigmaSst,   fileNameOutArcSigmaPod1,   fileNameOutArcSigmaPod2;
    FileName fileNameOutEpochSigmaSst, fileNameOutEpochSigmaPod1, fileNameOutEpochSigmaPod2;
    FileName fileNameOutCovSst,        fileNameOutCovPod1,        fileNameOutCovPod2;
    FileName fileNameOutSigmasCovSst;
    FileName fileNameResidualsSst,     fileNameResidualsPod1,     fileNameResidualsPod2;

    FileName fileNameArcList;
    FileName fileNameInArcSigmaSst,   fileNameInArcSigmaPod1,   fileNameInArcSigmaPod2;
    FileName fileNameInEpochSigmaSst, fileNameInEpochSigmaPod1, fileNameInEpochSigmaPod2;
    FileName fileNameInCovSst,        fileNameInCovPod1,        fileNameInCovPod2;
    FileName fileNameInSigmasCovSst;
    std::vector<FileName> fileNamesCovSst;

    FileName    fileNameInCovEpochPod1,  fileNameInCovEpochPod2;
    Double      sigma0Sst=1, sigma0Pod1=1, sigma0Pod2=1;
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
    if(readConfigSequence(config, "estimateSstArcCovarianceSigmas", Config::OPTIONAL, "", ""))
    {
      estimateSigmasCovSst = TRUE;
      readConfig(config, "outputfileSigmasCovarianceMatrixArc", fileNameOutSigmasCovSst, Config::OPTIONAL, "", "one variance factor per matrix");
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
      readConfig(config, "sigma",                              sigma0Sst,               Config::DEFAULT,  "1", "apriori factor of covariance function");
      readConfig(config, "inputfileSigmasPerArc",              fileNameInArcSigmaSst,   Config::OPTIONAL, "",  "apriori different accuaries for each arc (multiplicated with sigma)");
      readConfig(config, "inputfileSigmasPerEpoch",            fileNameInEpochSigmaSst, Config::OPTIONAL, "",  "apriori different accuaries for each epoch");
      readConfig(config, "inputfileCovarianceFunction",        fileNameInCovSst,        Config::OPTIONAL, "",  "approximate covariances in time");
      readConfig(config, "inputfileCovarianceMatrixArc",       fileNamesCovSst,         Config::OPTIONAL, "",  "Must be given per sst arc with correct dimensions.");
      readConfig(config, "inputfileSigmasCovarianceMatrixArc", fileNameInSigmasCovSst,  Config::OPTIONAL, "",  "Vector with one sigma for each <inputfileCovarianceMatrixArc>");
      readConfig(config, "sampling",                           samplingSst,             Config::DEFAULT,  "5", "[seconds] sampling of the covariance function");
      endSequence(config);
    }
    if(readConfigSequence(config, "covariancePod1", Config::MUSTSET, "", ""))
    {
      readConfig(config, "sigma",                        sigma0Pod1,               Config::DEFAULT,  "1",  "apriori factor of covariance function");
      readConfig(config, "inputfileSigmasPerArc",        fileNameInArcSigmaPod1,   Config::OPTIONAL, "",   "apriori different accuaries for each arc (multiplicated with sigma)");
      readConfig(config, "inputfileSigmasPerEpoch",      fileNameInEpochSigmaPod1, Config::OPTIONAL, "",   "apriori different accuaries for each epoch");
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
    readConfig(config, "downweightPod",          downweightPod,      Config::DEFAULT,  "1",    "downweight factor for POD");
    readConfig(config, "inputfileArcList",       fileNameArcList,    Config::OPTIONAL, "",     "list to correspond points of time to arc numbers");
    readConfig(config, "iterationCount",         iterCount,          Config::DEFAULT,  "3",    "(maximum) number of iterations for the estimation of calibration parameter and error PSD");
    readConfig(config, "variableNameIterations", iterVariableName,   Config::OPTIONAL, "",     "All output fileNames in preprocessing iteration are expanded with this variable prior to writing to disk");
    readConfig(config, "defaultBlockSize",       defaultBlockSize,   Config::DEFAULT,  "2048", "block size of static normal equation blocks");
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
    CovSst.resize(arcCount);
    if(fileNamesCovSst.size())
    {
      logStatus<<"read arc-wise sst covariance matrices"<<Log::endl;
      Parallel::forEachProcess(arcCount, [&](UInt arcNo)
      {
        CovSst.at(arcNo).resize(fileNamesCovSst.size());
        for(UInt i=0; i<fileNamesCovSst.size(); i++)
          readFileMatrix(fileNamesCovSst.at(i).appendBaseName(".arc"+arcNo%"%03i"s), CovSst.at(arcNo).at(i));
      }, processNo, comm);
    }

    sigmasCovSst = Vector(fileNamesCovSst.size(), 1.);
    if(!fileNameInSigmasCovSst.empty())
    {
      readFileMatrix(fileNameInSigmasCovSst, sigmasCovSst);
      if(sigmasCovSst.rows() != fileNamesCovSst.size())
        throw(Exception("Number of sigmas not compatible with number of given arc-wise SST covariance matrices"));
    }

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
    covFuncSst.column(1,1)  *= std::pow(sigma0Sst, 2);
    covFuncPod1.column(1,3) *= std::pow(sigma0Pod1,2);
    covFuncPod2.column(1,3) *= std::pow(sigma0Pod2,2);
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
          if((s2==s2) && (s2>0))
            sigma2ShortTimeModel = s2;
        }
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
        Parallel::forEachProcess(arcListResidualsSst,  [this](UInt arcNo) {return arcListResidualsSst.at(arcNo);},  processNo, comm);
        Parallel::forEachProcess(arcListResidualsPod1, [this](UInt arcNo) {return arcListResidualsPod1.at(arcNo);}, processNo, comm);
        Parallel::forEachProcess(arcListResidualsPod2, [this](UInt arcNo) {return arcListResidualsPod2.at(arcNo);}, processNo, comm);

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

      // compute redundancies
      // --------------------
      if((estimateArcSigmas || estimateCovarianceFunctionVCE || estimateSigmasCovSst))
      {
        logStatus<<"compute redundancies"<<Log::endl;
        sigmaSstNew = sigmaPod1New = sigmaPod2New = Vector(arcCount);
        ePeSst  = redundancySst  = Matrix(covLengthSst,  1);
        ePePod1 = redundancyPod1 = Matrix(covLengthPod1, 3); // for x,y,z
        ePePod2 = redundancyPod2 = Matrix(covLengthPod2, 3); // for x,y,z
        ePeCovSst = redundancyCovSst = Vector(sigmasCovSst.rows());
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
        Parallel::forEachProcess(arcListEpochSigmaSst,  [this](UInt arcNo) {return arcListEpochSigmaSst.at(arcNo);},  processNo, comm);
        Parallel::forEachProcess(arcListEpochSigmaPod1, [this](UInt arcNo) {return arcListEpochSigmaPod1.at(arcNo);}, processNo, comm);
        Parallel::forEachProcess(arcListEpochSigmaPod2, [this](UInt arcNo) {return arcListEpochSigmaPod2.at(arcNo);}, processNo, comm);

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
        } // if(Parallel::isMaster(comm))
        Parallel::broadCast(PsdSst,  0, comm);
        Parallel::broadCast(PsdPod1, 0, comm);
        Parallel::broadCast(PsdPod2, 0, comm);
        // compute new covariance function
        copy(CosTransformSst  * PsdSst,  covFuncSst.column(1,1));
        copy(CosTransformPod1 * PsdPod1, covFuncPod1.column(1,3));
        copy(CosTransformPod2 * PsdPod2, covFuncPod2.column(1,3));
      } // if(estimateCovarianceFunctions)

      // =============================================

      // estimate variance factor for arc-wise SST covariance matrices
      // -------------------------------------------------------------
      if(estimateSigmasCovSst)
      {
        Parallel::reduceSum(ePeCovSst, 0, comm);
        Parallel::reduceSum(redundancyCovSst, 0, comm);
        if(Parallel::isMaster(comm))
        {
          for(UInt i=0; i<sigmasCovSst.rows(); i++)
          {
            const Double alpha = std::sqrt(ePeCovSst(i)/redundancyCovSst(i));
            sigmasCovSst(i) *= alpha;
            logStatus<<"  sigma of SST covariance matrix (current/total): "<<alpha<<"/"<<sigmasCovSst(i)<<Log::endl;
          }
        }
        Parallel::broadCast(sigmasCovSst, 0, comm);
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
      if(Parallel::isMaster(comm) && !fileNameOutSigmasCovSst.empty())
      {
        logStatus<<"write arc-wise SST variance factors <"<<fileNameOutSigmasCovSst(variableIteration)<<">"<<Log::endl;
        writeFileMatrix(fileNameOutSigmasCovSst(variableIteration), sigmasCovSst);
      }

      PsdPod1                 *= std::pow(downweightPod, 2);
      PsdPod2                 *= std::pow(downweightPod, 2);
      covFuncPod1.column(1,3) *= std::pow(downweightPod, 2);
      covFuncPod2.column(1,3) *= std::pow(downweightPod, 2);
    } // for(iter)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

UInt PreprocessingSst::findInterval(UInt arcNo) const
{
  for(UInt idInterval=0; idInterval+1<arcsInterval.size(); idInterval++)
    if(arcNo < arcsInterval.at(idInterval+1))
      return idInterval;
  return 0;
}

/***********************************************/

void PreprocessingSst::decorrelate(UInt arcNo, UInt countSst, UInt countPod1, UInt countPod2,
                                   Matrix &WSst, Matrix &WPod1, Matrix &WPod2, const std::list<MatrixSlice> &A)
{
  try
  {
    UInt obsCount = 0;
    const UInt idxSst  = obsCount; obsCount += countSst;
    const UInt idxPod1 = obsCount; obsCount += 3*countPod1;
    const UInt idxPod2 = obsCount; obsCount += 3*countPod2;

    if(countSst)
    {
      std::list<MatrixSlice> WA;
      for(const MatrixSlice &a : A)
        WA.push_back(a.row(idxSst, countSst));

      WSst = Matrix();
      if(CovSst.at(arcNo).size())
      {
        WSst = std::pow(sigmasCovSst.at(0),2) * CovSst.at(arcNo).at(0);
        for(UInt i=1; i<CovSst.at(arcNo).size(); i++)
          axpy(std::pow(sigmasCovSst.at(i),2), CovSst.at(arcNo).at(i), WSst);
      }

      CovarianceSst::decorrelate(observationArc.at(arcNo).timesSst, sigmaSst(arcNo), arcListEpochSigmaSst.at(arcNo), covFuncSst, WSst, WA);
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

void PreprocessingSst::buildNormals(UInt arcNo)
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
    Matrix WSst, WPod1, WPod2;
    decorrelate(arcNo, observationArc.at(arcNo).timesSst.size(),
                observationArc.at(arcNo).timesPod1.size(), observationArc.at(arcNo).timesPod2.size(),
                WSst, WPod1, WPod2, {Wl, WA, WB});

    normals.accumulate(findInterval(arcNo), Wl, WA, WB);
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
    Matrix Wl = observationArc.at(arcNo).l;
    Matrix WA = observationArc.at(arcNo).A;
    Matrix WB = observationArc.at(arcNo).B;
    Matrix WSst, WPod1, WPod2;
    decorrelate(arcNo, countSst, countPod1, countPod2, WSst, WPod1, WPod2, {Wl, WA, WB});

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
    normals.designMatMult(findInterval(arcNo), -1., WA, x,  We);
    normals.designMatMult(findInterval(arcNo), +1., WA, Wz, WAz);

    // ============================================

    if(!estimateCovarianceFunctionVCE)
    {
      if(countSst)
      {
        const Double redundancy = countSst - quadsum(WAz.row(idxSst, countSst)) - quadsum(WB.row(idxSst, countSst));
        sigmaSstNew(arcNo) = std::sqrt(quadsum(We.row(idxSst,  countSst))/redundancy) * sigmaSst(arcNo);
      }
      if(countPod1)
      {
        const Double redundancy = 3*countPod1 - quadsum(WAz.row(idxPod1, 3*countPod1)) - quadsum(WB.row(idxPod1, 3*countPod1));
        sigmaPod1New(arcNo) = std::sqrt(quadsum(We.row(idxPod1,  3*countPod1))/redundancy) * sigmaPod1(arcNo);
      }
      if(countPod2)
      {
        const Double redundancy = 3*countPod2 - quadsum(WAz.row(idxPod2, 3*countPod2)) - quadsum(WB.row(idxPod2, 3*countPod2));
        sigmaPod2New(arcNo) = std::sqrt(quadsum(We.row(idxPod2,  3*countPod2))/redundancy) * sigmaPod2(arcNo);
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
      if(estimateSigmasCovSst)
        for(UInt i=0; i<CovSst.at(arcNo).size(); i++)
          Vce::matrix(R, WWe, std::pow(sigmasCovSst(i),2) * CovSst.at(arcNo).at(i), ePeCovSst(i), redundancyCovSst(i));
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
    Matrix e = observationArc.at(arcNo).l;
    normals.designMatMult(findInterval(arcNo), -1., observationArc.at(arcNo).A, x, e);

    // eliminate arc dependent parameters
    // ----------------------------------
    if(observationArc.at(arcNo).B.size())
    {
      Matrix We = e;
      Matrix WB = observationArc.at(arcNo).B;
      Matrix WSst, WPod1, WPod2;
      decorrelate(arcNo, countSst, countPod1, countPod2, WSst, WPod1, WPod2, {We, WB});
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
      const Double e2 = std::pow(arcListResidualsSst.at(arcNo).at(i).rangeRate, 2);
      arcListEpochSigmaSst.at(arcNo).at(i).sigma = 0.;
      if(e2 > thresholdSst)
        arcListEpochSigmaSst.at(arcNo).at(i).sigma = std::sqrt(e2-thresholdSst);
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
/***********************************************/
