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

  static constexpr UInt TYPESIZE = 3;
  enum Type : UInt {SST, POD1, POD2};
  std::array<std::string, TYPESIZE> typeName; // = {"SST", "POD1", "POD2"};

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
  std::array<UInt,   TYPESIZE> covColumns;
  std::array<Vector, TYPESIZE> sigma, sigmaNew;
  std::array<Matrix, TYPESIZE> covFunc, Psd, ePe, redundancy, CosTransform;
  std::array<Double, TYPESIZE> sampling;
  std::array<std::vector<ObservationSigmaArc>, TYPESIZE> arcListEpochSigma;

  std::vector<std::vector<Matrix>> CovSst; // Several independant matrices per arc
  Vector  sigmasCovSst, ePeCovSst, redundancyCovSst;

  // residuals
  std::array<std::vector<Arc>, TYPESIZE> arcListResiduals;

  UInt findInterval(UInt arcNo) const;
  void decorrelate(UInt arcNo, const std::array<UInt, TYPESIZE> &count, std::array<Matrix, TYPESIZE> &W, const std::list<MatrixSlice> &A);
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
    typeName   = {"SST", "POD1", "POD2"};
    covColumns = {1, 3, 3};
    estimateCovarianceFunctionVCE = estimateArcSigmas = estimateEpochSigmas = estimateResiduals = FALSE;
    estimateSigmaShortTimeModel = FALSE;
    estimateSigmasCovSst = FALSE;

    FileName                       fileNameSolution, fileNameSigmax, fileNameParaName;
    std::array<FileName, TYPESIZE> fileNameOutArcSigma, fileNameOutEpochSigma, fileNameOutCovFunc, fileNameResiduals;
    std::array<FileName, TYPESIZE> fileNameInArcSigma, fileNameInEpochSigma, fileNameInCovFunc;
    FileName                       fileNameOutSigmasCovSst, fileNameInSigmasCovSst;
    FileName                       fileNameInCovEpochPod1,  fileNameInCovEpochPod2;
    std::vector<FileName>          fileNamesCovSst;
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
      readConfig(config, "outputfileSigmasPerArcSst",  fileNameOutArcSigma.at(SST),  Config::OPTIONAL, "", "accuracies of each arc (SST)");
      readConfig(config, "outputfileSigmasPerArcPod1", fileNameOutArcSigma.at(POD1), Config::OPTIONAL, "", "accuracies of each arc (POD1)");
      readConfig(config, "outputfileSigmasPerArcPod2", fileNameOutArcSigma.at(POD2), Config::OPTIONAL, "", "accuracies of each arc (POD2)");
      endSequence(config);
    }
    if(readConfigSequence(config, "estimateEpochSigmas", Config::OPTIONAL, "", ""))
    {
      estimateEpochSigmas = TRUE;
      readConfig(config, "outputfileSigmasPerEpochSst",  fileNameOutEpochSigma.at(SST),  Config::OPTIONAL, "", "accuracies of each epoch (SST)");
      readConfig(config, "outputfileSigmasPerEpochPod1", fileNameOutEpochSigma.at(POD1), Config::OPTIONAL, "", "accuracies of each epoch (POD1)");
      readConfig(config, "outputfileSigmasPerEpochPod2", fileNameOutEpochSigma.at(POD2), Config::OPTIONAL, "", "accuracies of each epoch (POD2)");
      endSequence(config);
    }
    if(readConfigSequence(config, "estimateCovarianceFunctions", Config::OPTIONAL, "", ""))
    {
      estimateCovarianceFunctionVCE = TRUE;
      readConfig(config, "outputfileCovarianceFunctionSst",  fileNameOutCovFunc.at(SST),  Config::OPTIONAL, "", "covariance function");
      readConfig(config, "outputfileCovarianceFunctionPod1", fileNameOutCovFunc.at(POD1), Config::OPTIONAL, "", "covariance functions for along, cross, radial direction");
      readConfig(config, "outputfileCovarianceFunctionPod2", fileNameOutCovFunc.at(POD2), Config::OPTIONAL, "", "covariance functions for along, cross, radial direction");
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
      readConfig(config, "outputfileSstResiduals",  fileNameResiduals.at(SST),  Config::OPTIONAL, "", "");
      readConfig(config, "outputfilePod1Residuals", fileNameResiduals.at(POD1), Config::OPTIONAL, "", "");
      readConfig(config, "outputfilePod2Residuals", fileNameResiduals.at(POD2), Config::OPTIONAL, "", "");
      endSequence(config);
    }
    readConfig(config, "observation", observationMisc, Config::MUSTSET,  "", "");
    if(readConfigSequence(config, "covarianceSst", Config::MUSTSET, "", ""))
    {
      readConfig(config, "sigma",                              sigma0.at(SST),               Config::DEFAULT,  "1", "apriori factor of covariance function");
      readConfig(config, "inputfileSigmasPerArc",              fileNameInArcSigma.at(SST),   Config::OPTIONAL, "",  "apriori different accuaries for each arc (multiplicated with sigma)");
      readConfig(config, "inputfileSigmasPerEpoch",            fileNameInEpochSigma.at(SST), Config::OPTIONAL, "",  "apriori different accuaries for each epoch");
      readConfig(config, "inputfileCovarianceFunction",        fileNameInCovFunc.at(SST),    Config::OPTIONAL, "",  "approximate covariances in time");
      readConfig(config, "inputfileCovarianceMatrixArc",       fileNamesCovSst,              Config::OPTIONAL, "",  "Must be given per sst arc with correct dimensions.");
      readConfig(config, "inputfileSigmasCovarianceMatrixArc", fileNameInSigmasCovSst,       Config::OPTIONAL, "",  "Vector with one sigma for each <inputfileCovarianceMatrixArc>");
      readConfig(config, "sampling",                           sampling.at(SST),             Config::DEFAULT,  "5", "[seconds] sampling of the covariance function");
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
    for(UInt idType : {SST, POD1, POD2})
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

    std::array<UInt, TYPESIZE> covLength; covLength.fill(0);
    for(UInt idType : {SST, POD1, POD2})
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
      covFunc.at(idType)  = Vce::readCovarianceFunction(fileNameInCovFunc.at(idType),  covLength.at(idType),  covColumns.at(idType), sampling.at(idType));
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
        for(UInt idType : {SST, POD1, POD2})
        {
          arcListResiduals.at(idType).clear();
          arcListResiduals.at(idType).resize(arcCount);
        }
        Parallel::forEachProcess(arcCount, [this](UInt arcNo) {computeResiduals(arcNo);}, processNo, comm, FALSE/*timing*/);
        for(UInt idType : {SST, POD1, POD2})
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
      if((estimateArcSigmas || estimateCovarianceFunctionVCE || estimateSigmasCovSst))
      {
        logStatus<<"compute redundancies"<<Log::endl;
        for(UInt idType : {SST, POD1, POD2})
        {
          sigmaNew.at(idType) = Vector(arcCount);
          ePe.at(idType) = redundancy.at(idType) = Matrix(covLength.at(idType), covColumns.at(idType));
        }
        ePeCovSst = redundancyCovSst = Vector(sigmasCovSst.rows());
        Parallel::forEachProcess(arcCount, [this](UInt arcNo) {computeRedundancies(arcNo);}, processNo, comm);
      }

      // sigmas per arc
      // --------------
      if(estimateArcSigmas)
      {
        for(UInt idType : {SST, POD1, POD2})
        {
          Parallel::reduceSum(sigmaNew.at(idType), 0, comm);
          if(Parallel::isMaster(comm))
          {
            sigma.at(idType) = (1./Vce::meanSigma(sigmaNew.at(idType))) * sigmaNew.at(idType);
            logInfo<<"  "<<typeName.at(idType)<<" sigma per arc (median): "<<Vce::meanSigma(sigmaNew.at(idType))<<Log::endl;
          }
          Parallel::broadCast(sigma.at(idType), 0, comm);
        }

        for(UInt idType : {SST, POD1, POD2})
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
        for(UInt idType : {SST, POD1, POD2})
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
        for(UInt idType : {SST, POD1, POD2})
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
        for(UInt idType : {SST, POD1, POD2})
          if(Parallel::isMaster(comm) && !fileNameOutCovFunc.at(idType).empty())
          {
            logStatus<<"write covariance function file <"<<fileNameOutCovFunc.at(idType)(variableIteration)<<">"<<Log::endl;
            writeFileMatrix(fileNameOutCovFunc.at(idType)(variableIteration), covFunc.at(idType));
          }
      }

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

        if(Parallel::isMaster(comm) && !fileNameOutSigmasCovSst.empty())
        {
          logStatus<<"write arc-wise SST variance factors <"<<fileNameOutSigmasCovSst(variableIteration)<<">"<<Log::endl;
          writeFileMatrix(fileNameOutSigmasCovSst(variableIteration), sigmasCovSst);
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

UInt PreprocessingSst::findInterval(UInt arcNo) const
{
  for(UInt idInterval=0; idInterval+1<arcsInterval.size(); idInterval++)
    if(arcNo < arcsInterval.at(idInterval+1))
      return idInterval;
  return 0;
}

/***********************************************/

void PreprocessingSst::decorrelate(UInt arcNo, const std::array<UInt, TYPESIZE> &count, std::array<Matrix, TYPESIZE> &W, const std::list<MatrixSlice> &A)
{
  try
  {
    // count observations and calculate index
    // --------------------------------------
    std::array<UInt, TYPESIZE> idx;
    std::array<std::list<MatrixSlice>, TYPESIZE> WA; // type specific slices of A
    UInt obsCount = 0;
    for(UInt idType : {SST, POD1, POD2})
    {
      idx.at(idType) = obsCount;
      obsCount      += covColumns.at(idType) * count.at(idType);
      for(const MatrixSlice &a : A)
        WA.at(idType).push_back(a.row(idx.at(idType), covColumns.at(idType) * count.at(idType)));
    }

    if(count.at(SST))
    {
      W.at(SST) = Matrix();
      if(CovSst.at(arcNo).size())
      {
        W.at(SST) = std::pow(sigmasCovSst.at(0),2) * CovSst.at(arcNo).at(0);
        for(UInt i=1; i<CovSst.at(arcNo).size(); i++)
          axpy(std::pow(sigmasCovSst.at(i),2), CovSst.at(arcNo).at(i), W.at(SST));
      }

      CovarianceSst::decorrelate(observationArc.at(arcNo).times.at(SST), sigma.at(SST)(arcNo), arcListEpochSigma.at(SST).at(arcNo), covFunc.at(SST), W.at(SST), WA.at(SST));
    }

    if(count.at(POD1))
      W.at(POD1) = CovariancePod::decorrelate(observationArc.at(arcNo).pod1, sigma.at(POD1)(arcNo), arcListEpochSigma.at(POD1).at(arcNo), fileCovEpochPod1.readArc(arcNo), covFunc.at(POD1), WA.at(POD1));

    if(count.at(POD2))
      W.at(POD2) = CovariancePod::decorrelate(observationArc.at(arcNo).pod2, sigma.at(POD2)(arcNo), arcListEpochSigma.at(POD2).at(arcNo), fileCovEpochPod2.readArc(arcNo), covFunc.at(POD2), WA.at(POD2));
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
    std::array<Matrix, TYPESIZE> W;
    decorrelate(arcNo, {observationArc.at(arcNo).times.at(SST).size(), observationArc.at(arcNo).times.at(POD1).size(), observationArc.at(arcNo).times.at(POD2).size()},
                W, {Wl, WA, WB});

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
    std::array<UInt, TYPESIZE> count, idx;
    UInt obsCount = 0;
    for(UInt idType : {SST, POD1, POD2})
    {
      count.at(idType) = observationArc.at(arcNo).times.at(idType).size();
      idx.at(idType)   = obsCount;
      obsCount        += covColumns.at(idType) * count.at(idType);
    }

    // Decorrelation
    // -------------
    Matrix Wl = observationArc.at(arcNo).l;
    Matrix WA = observationArc.at(arcNo).A;
    Matrix WB = observationArc.at(arcNo).B;
    std::array<Matrix, TYPESIZE> W;
    decorrelate(arcNo, count, W, {Wl, WA, WB});

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
      for(UInt idType : {SST, POD1, POD2})
        if(count.at(idType))
        {
          const Double redundancy = covColumns.at(idType)*count.at(idType)
                                  - quadsum(WAz.row(idx.at(idType), covColumns.at(idType)*count.at(idType)))
                                  - quadsum(WB.row (idx.at(idType), covColumns.at(idType)*count.at(idType)));
          if(redundancy > 0.3)
            sigmaNew.at(idType)(arcNo) = std::sqrt(quadsum(We.row(idx.at(idType),  covColumns.at(idType)*count.at(idType)))/redundancy) * sigma.at(idType)(arcNo);
        }
      return;
    } // if(!estimateCovarianceFunctionVCE)

    // ============================================

    // variance component estimation
    // -----------------------------
    Matrix R;
    Vector WWe;
    std::array<std::vector<UInt>, TYPESIZE> index;
    std::array<Double, TYPESIZE> ePeSum, redundancySum; ePeSum.fill(0); redundancySum.fill(0);
    for(UInt idType : {SST, POD1, POD2})
    {
      index.at(idType).resize(observationArc.at(arcNo).times.at(idType).size());
      for(UInt i=0; i<index.at(idType).size(); i++)
        index.at(idType).at(i) = static_cast<UInt>(std::round((observationArc.at(arcNo).times.at(idType).at(i)-observationArc.at(arcNo).times.at(idType).front()).seconds()/sampling.at(idType)));
    }

    if(count.at(SST))
    {
      Vce::redundancy(W.at(SST), We.row(idx.at(SST), count.at(SST)),
                      WAz.row(idx.at(SST), count.at(SST)), WB.row(idx.at(SST), count.at(SST)), R, WWe);
      Vce::psd(R, WWe, index.at(SST), sigma.at(SST)(arcNo), CosTransform.at(SST), Psd.at(SST),
               ePe.at(SST), redundancy.at(SST), ePeSum.at(SST), redundancySum.at(SST));

      // Arc-wise covariance matrices
      if(estimateSigmasCovSst)
        for(UInt i=0; i<CovSst.at(arcNo).size(); i++)
          Vce::matrix(R, WWe, std::pow(sigmasCovSst(i),2) * CovSst.at(arcNo).at(i), ePeCovSst(i), redundancyCovSst(i));
    }

    for(UInt idType : {POD1, POD2})
      if(count.at(idType))
      {
        Vce::redundancy(W.at(idType), We.row(idx.at(idType), 3*count.at(idType)),
                        WAz.row(idx.at(idType), 3*count.at(idType)), WB.row(idx.at(idType), 3*count.at(idType)), R, WWe);
        Vce::psd(R, WWe, index.at(idType), sigma.at(idType)(arcNo), CosTransform.at(idType), Psd.at(idType),
                 ePe.at(idType), redundancy.at(idType), ePeSum.at(idType), redundancySum.at(idType));
      }


    for(UInt idType : {SST, POD1, POD2})  // compute new sigma (for this arc)
      if(redundancySum.at(idType) > 0.3)
        sigmaNew.at(idType)(arcNo) = std::sqrt(ePeSum.at(idType)/redundancySum.at(idType)) * sigma.at(idType)(arcNo);
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
    std::array<UInt, TYPESIZE> count, idx;
    UInt obsCount = 0;
    for(UInt idType : {SST, POD1, POD2})
    {
      count.at(idType) = observationArc.at(arcNo).times.at(idType).size();
      idx.at(idType)   = obsCount;
      obsCount        += covColumns.at(idType) * count.at(idType);
    }

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
      std::array<Matrix, TYPESIZE> W;
      decorrelate(arcNo, count, W, {We, WB});
      Vector tau = QR_decomposition(WB);
      QTransMult(WB, tau, We); // transform observations: l:= Q'l
      Matrix y = We.row(0, tau.rows());
      triangularSolve(1., WB.row(0, tau.rows()), y);
      matMult(-1, observationArc.at(arcNo).B, y, e);
    }

    // create Sst arc
    // --------------
    for(UInt i=0; i<count.at(SST); i++)
    {
      SatelliteTrackingEpoch epoch;
      epoch.time  = observationArc.at(arcNo).times.at(SST).at(i);
      epoch.range = epoch.rangeRate = epoch.rangeAcceleration = e(idx.at(SST)+i,0);
      arcListResiduals.at(SST).at(arcNo).push_back(epoch);
    }

    // create Pod arcs
    // ---------------
    for(UInt idType : {POD1, POD2})
      for(UInt i=0; i<count.at(idType); i++)
      {
        OrbitEpoch epoch;
        epoch.time     = observationArc.at(arcNo).times.at(idType).at(i);
        epoch.position = Vector3d(e.row(idx.at(idType)+3*i, 3));
        arcListResiduals.at(idType).at(arcNo).push_back(epoch);
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

    if(covFunc.at(SST).size())
    {
      const Double threshold2 = std::pow(huber*sigma.at(SST)(arcNo), 2) * covFunc.at(SST)(0,1);
      for(UInt i=0; i<observationArc.at(arcNo).times.at(SST).size(); i++)
      {
        const Double e2 = std::pow(dynamic_cast<const SatelliteTrackingEpoch &>(arcListResiduals.at(SST).at(arcNo).at(i)).rangeRate, 2);
        arcListEpochSigma.at(SST).at(arcNo).at(i).sigma = 0.;
        if(e2 > threshold2)
          arcListEpochSigma.at(SST).at(arcNo).at(i).sigma = std::sqrt(e2-threshold2);
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
