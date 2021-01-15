/***********************************************/
/**
* @file preprocessingVariationalEquationSstFit.cpp
*
* @brief Fit variational equations to sst/orbit observations.
*
* @author Torsten Mayer-Guerr
* @date 2015-07-06
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program fits two \configFile{inputfileVariational1/2}{variationalEquation} to satellite-to-satellite-tracking (SST) and orbit
observations in a GRACE like configuration. It works similar to \program{PreprocessingVariationalEquationOrbitFit},
see there for details.

As the relative weighting of the observation types is important complex description of the covariances can be set with
\configClass{covarianceSst}{covarianceSstType}, \configClass{covariancePod1}{covariancePodType}, \configClass{covariancePod2}{covariancePodType}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"
#include "files/fileVariationalEquation.h"
#include "files/fileParameterName.h"
#include "misc/varianceComponentEstimation.h"
#include "misc/observation/observationMiscSstVariational.h"
#include "misc/observation/covariancePod.h"
#include "misc/observation/covarianceSst.h"
#include "misc/observation/variationalEquation.h"
#include "misc/observation/variationalEquationFromFile.h"

/***** CLASS ***********************************/

/** @brief Fit variational equations to sst/orbit observations.
* @ingroup programsGroup */
class PreprocessingVariationalEquationSstFit
{
public:
  ObservationMiscSstVariationalPtr observationMisc;
  CovarianceSstPtr covSst;
  CovariancePodPtr covPod1, covPod2;
  std::vector<ObservationMiscSst::Arc> observationArc;

  // normal equations
  // ----------------
  Matrix N_Sst, N_Pod;              // =A'PA, Normal matrix
  Vector n_Sst, n_Pod;              // =A'Pl, right hand side
  Double lPl_Sst, lPl_Pod;          // =l'Pl, weighted norm of the observations
  UInt   obsCountSst, obsCountPod;  // number of observations
  UInt   outlierCountSst, outlierCountPod;
  Double sigma0Sst, sigma0Pod;
  Vector x, x1, x2;

  void run(Config &config, Parallel::CommunicatorPtr comm);

  void buildNormals(UInt arcNo);
};

GROOPS_REGISTER_PROGRAM(PreprocessingVariationalEquationSstFit, PARALLEL, "fit variational equations to sst/orbit observations", Preprocessing, VariationalEquation)

/***********************************************/

void PreprocessingVariationalEquationSstFit::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName fileNameOutVariational1, fileNameOutVariational2;
    FileName fileNameOutOrbit1, fileNameOutOrbit2;
    FileName fileNameOutSolution1, fileNameOutSolution2;
    UInt     iterCount;

    readConfig(config, "outputfileVariational1", fileNameOutVariational1, Config::MUSTSET,  "", "approximate position and integrated state matrix");
    readConfig(config, "outputfileVariational2", fileNameOutVariational2, Config::MUSTSET,  "", "approximate position and integrated state matrix");
    readConfig(config, "outputfileOrbit1",       fileNameOutOrbit1,       Config::OPTIONAL, "", "integrated orbit");
    readConfig(config, "outputfileOrbit2",       fileNameOutOrbit2,       Config::OPTIONAL, "", "integrated orbit");
    readConfig(config, "outputfileSolution1",    fileNameOutSolution1,    Config::OPTIONAL, "", "estimated calibration and state parameters");
    readConfig(config, "outputfileSolution2",    fileNameOutSolution2,    Config::OPTIONAL, "", "estimated calibration and state parameters");
    observationMisc = ObservationMiscSstVariationalPtr(new ObservationMiscSstVariational(config));
    readConfig(config, "covarianceSst",          covSst,                  Config::MUSTSET,  "", "covariance matrix of satellite to satellite tracking observations");
    readConfig(config, "covariancePod1",         covPod1,                 Config::MUSTSET,  "", "covariance matrix of kinematic orbits (satellite 1)");
    readConfig(config, "covariancePod2",         covPod2,                 Config::MUSTSET,  "", "covariance matrix of kinematic orbits (satellite 2)");
    readConfig(config, "iterationCount",         iterCount,               Config::DEFAULT,  "10", "for the estimation of calibration parameter and error PSD");
    if(isCreateSchema(config)) return;

    // =============================================

    logStatus<<"set up observation equations"<<Log::endl;
    observationArc.resize(observationMisc->arcCount());
    std::vector<UInt> processNo = Parallel::forEach(observationMisc->arcCount(), [this](UInt arcNo) {observationArc.at(arcNo) = observationMisc->computeArc(arcNo, covSst, covPod1, covPod2);}, comm);

    x = Vector(observationMisc->parameterCount());
    sigma0Sst = 1;
    sigma0Pod = 1;
    for(UInt iter=0; iter<iterCount; iter++)
    {
      // build normals
      // -------------
      logStatus<<"accumulate normal equations"<<Log::endl;
      logInfo<<"  parameter count = "<<observationMisc->parameterCount()<<Log::endl;
      N_Sst   = Matrix(observationMisc->parameterCount(), Matrix::SYMMETRIC);
      N_Pod   = Matrix(observationMisc->parameterCount(), Matrix::SYMMETRIC);
      n_Sst   = Vector(observationMisc->parameterCount());
      n_Pod   = Vector(observationMisc->parameterCount());
      lPl_Sst = 0;
      lPl_Pod = 0;
      obsCountSst = outlierCountSst = 0;
      obsCountPod = outlierCountPod = 0;

      Parallel::forEachProcess(observationMisc->arcCount(), [this](UInt arcNo) {buildNormals(arcNo);}, processNo, comm);

      Parallel::reduceSum(N_Sst,           0, comm);
      Parallel::reduceSum(N_Pod,           0, comm);
      Parallel::reduceSum(n_Sst,           0, comm);
      Parallel::reduceSum(n_Pod,           0, comm);
      Parallel::reduceSum(lPl_Sst,         0, comm);
      Parallel::reduceSum(lPl_Pod,         0, comm);
      Parallel::reduceSum(obsCountSst,     0, comm);
      Parallel::reduceSum(obsCountPod,     0, comm);
      Parallel::reduceSum(outlierCountSst, 0, comm);
      Parallel::reduceSum(outlierCountPod, 0, comm);
      UInt outlierCount = outlierCountSst + outlierCountPod;

      // Estimate parameters
      // -------------------
      if(Parallel::isMaster(comm))
      {
        Matrix N = pow(1./sigma0Sst,2) * N_Sst + pow(1./sigma0Pod,2) * N_Pod;
        Vector n = pow(1./sigma0Sst,2) * n_Sst + pow(1./sigma0Pod,2) * n_Pod;

        // Regularize unused parameters
        for(UInt i=0; i<N.rows(); i++)
          if(N(i,i) == 0)
            N(i,i) = 1.0;

        logStatus<<"solve system of normal equations"<<Log::endl;
        x = solve(N,n);

        // variance component estimation
        cholesky2Inverse(N);
        fillSymmetric(N);
        fillSymmetric(N_Pod);
        fillSymmetric(N_Sst);
        const Double ePe_Sst = lPl_Sst - 2*inner(x, n_Sst) + inner(x, N_Sst*x);
        const Double ePe_Pod = lPl_Pod - 2*inner(x, n_Pod) + inner(x, N_Pod*x);
        const Double r_Sst   = obsCountSst - pow(1./sigma0Sst,2)*inner(N, N_Sst);
        const Double r_Pod   = obsCountPod - pow(1./sigma0Pod,2)*inner(N, N_Pod);
        sigma0Sst = Vce::standardDeviation(ePe_Sst, r_Sst, 2.5/*huber*/, 1./*huberPower*/);
        sigma0Pod = Vce::standardDeviation(ePe_Pod, r_Pod, 2.5/*huber*/, 1./*huberPower*/);

        logInfo<<"  aposteriori sigma (SST): "<<sigma0Sst<<Log::endl;
        logInfo<<"  aposteriori sigma (POD): "<<sigma0Pod<<Log::endl;
        logInfo<<"  outlier (SST) "<<outlierCountSst<<" of "<<obsCountSst<<" ("<<100.*outlierCountSst/obsCountSst<<"%)"<<Log::endl;
        logInfo<<"  outlier (POD) "<<outlierCountPod<<" of "<<obsCountPod<<" ("<<100.*outlierCountPod/obsCountPod<<"%)"<<Log::endl;
      }

      Parallel::broadCast(outlierCount, 0, comm);
      if((iter>0) && (outlierCount==0))
        break;

      Parallel::broadCast(x, 0, comm);
      Parallel::broadCast(sigma0Sst, 0, comm);
      Parallel::broadCast(sigma0Pod, 0, comm);
    } // for(iter)

    // Split solution
    // --------------
    const UInt gravityCount = observationMisc->variationalEquation1.parameterCountGravity();
    const UInt state1Count  = observationMisc->variationalEquation1.parameterCount() - gravityCount;
    const UInt state2Count  = observationMisc->variationalEquation2.parameterCount() - gravityCount;

    UInt countAParameter = 0;
    const UInt idxGravity  = countAParameter; countAParameter += gravityCount;
    const UInt idxState1   = countAParameter; countAParameter += state1Count;
    const UInt idxState2   = countAParameter; countAParameter += state2Count;

    x1 = Vector(gravityCount+state1Count);
    if(gravityCount) copy(x.row(idxGravity, gravityCount), x1.row(0, gravityCount));
    if(state1Count)  copy(x.row(idxState1, state1Count),   x1.row(gravityCount, state1Count));

    x2 = Vector(gravityCount+state2Count);
    if(gravityCount) copy(x.row(idxGravity, gravityCount), x2.row(0, gravityCount));
    if(state2Count)  copy(x.row(idxState2, state2Count),   x2.row(gravityCount, state2Count));

    if(Parallel::isMaster(comm) && !fileNameOutSolution1.empty())
    {
      logStatus<<"write solution to <"<<fileNameOutSolution1<<">"<<Log::endl;
      writeFileMatrix(fileNameOutSolution1, x1);

      std::vector<ParameterName> parameterName;
      observationMisc->variationalEquation1.parameterName(parameterName);
      writeFileParameterName(fileNameOutSolution1.replaceFullExtension(".parameterName.txt"), parameterName);
    }

    if(Parallel::isMaster(comm) && !fileNameOutSolution2.empty())
    {
      logStatus<<"write solution to <"<<fileNameOutSolution2<<">"<<Log::endl;
      writeFileMatrix(fileNameOutSolution2, x2);

      std::vector<ParameterName> parameterName;
      observationMisc->variationalEquation2.parameterName(parameterName);
      writeFileParameterName(fileNameOutSolution2.replaceFullExtension(".parameterName.txt"), parameterName);
    }

    // =============================================================================

    // Improve orbits with estimated parameters
    // ----------------------------------------
    logStatus<<"Improve orbits with estimated parameters"<<Log::endl;
    std::vector<VariationalEquationArc> arcs1(observationMisc->variationalEquation1.arcCount());
    std::vector<VariationalEquationArc> arcs2(observationMisc->variationalEquation2.arcCount());
    Parallel::forEach(arcs1, [&](UInt arcNo) {return observationMisc->variationalEquation1.refineVariationalEquationArc(arcNo, x1);}, comm);
    Parallel::forEach(arcs2, [&](UInt arcNo) {return observationMisc->variationalEquation2.refineVariationalEquationArc(arcNo, x2);}, comm);

    // =============================================================================

    if(Parallel::isMaster(comm) && !fileNameOutVariational1.empty())
    {
      logStatus<<"write variational equation to file <"<<fileNameOutVariational1<<">"<<Log::endl;
      writeFileVariationalEquation(fileNameOutVariational2, observationMisc->variationalEquation1.satellite(), arcs1);
    }

    if(Parallel::isMaster(comm) && !fileNameOutVariational2.empty())
    {
      logStatus<<"write variational equation to file <"<<fileNameOutVariational2<<">"<<Log::endl;
      writeFileVariationalEquation(fileNameOutVariational2, observationMisc->variationalEquation2.satellite(), arcs2);
    }

    if(Parallel::isMaster(comm) && !fileNameOutOrbit1.empty())
    {
      logStatus<<"write orbit to file <"<<fileNameOutOrbit1<<">"<<Log::endl;
      std::list<Arc> arcList;
      for(UInt arcNo=0; arcNo<arcs1.size(); arcNo++)
        arcList.push_back( arcs1.at(arcNo).orbitArc() );
      InstrumentFile::write(fileNameOutOrbit1, arcList);
    }

    if(Parallel::isMaster(comm) && !fileNameOutOrbit2.empty())
    {
      logStatus<<"write orbit to file <"<<fileNameOutOrbit2<<">"<<Log::endl;
      std::list<Arc> arcList;
      for(UInt arcNo=0; arcNo<arcs2.size(); arcNo++)
        arcList.push_back( arcs2.at(arcNo).orbitArc() );
      InstrumentFile::write(fileNameOutOrbit2, arcList);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void PreprocessingVariationalEquationSstFit::buildNormals(UInt arcNo)
{
  try
  {
    Vector l = observationArc.at(arcNo).l;
    Matrix A = observationArc.at(arcNo).A;

    // count observations and calculate index
    // --------------------------------------
    const UInt countSst  = observationArc.at(arcNo).timesSst.size();
    const UInt countPod1 = observationArc.at(arcNo).timesPod1.size();
    const UInt countPod2 = observationArc.at(arcNo).timesPod2.size();

    UInt obsCount = 0;
    const UInt idxSst  = obsCount; obsCount += countSst;
    const UInt idxPod1 = obsCount; obsCount += 3*countPod1;
    const UInt idxPod2 = obsCount; obsCount += 3*countPod2;

    // downweight outliers
    if(quadsum(x))
    {
      const Double huber = 2.5;

      Vector e = l;
      matMult(-1, A, x, e);

      for(UInt k=0; k<countSst; k++)
      {
        const Double s = fabs(e(idxSst+k));
        if(s>huber*sigma0Sst)
        {
          l.row(idxSst+k) *= huber*sigma0Sst/s;
          A.row(idxSst+k) *= huber*sigma0Sst/s;
          outlierCountSst++;
        }
      }

      for(UInt k=0; k<countPod1; k++)
      {
        const Double s = sqrt(quadsum(e.row(idxPod1+3*k,3))/3);
        if(s>huber*sigma0Pod)
        {
          l.row(idxPod1+3*k,3) *= huber*sigma0Pod/s;
          A.row(idxPod1+3*k,3) *= huber*sigma0Pod/s;
          outlierCountPod += 3;
        }
      }

      for(UInt k=0; k<countPod2; k++)
      {
        const Double s = sqrt(quadsum(e.row(idxPod2+3*k,3))/3);
        if(s>huber*sigma0Pod)
        {
          l.row(idxPod2+3*k,3) *= huber*sigma0Pod/s;
          A.row(idxPod2+3*k,3) *= huber*sigma0Pod/s;
          outlierCountPod += 3;
        }
      }

    } // if(quadsum(x))

    // Sst
    rankKUpdate(1., A.row(idxSst, countSst), N_Sst);
    matMult(1., A.row(idxSst, countSst).trans(), l.row(idxSst, countSst), n_Sst);
    lPl_Sst      += quadsum(l.row(idxSst, countSst));
    obsCountSst  += countSst;

    // Pod
    rankKUpdate(1., A.row(idxPod1, 3*(countPod1+countPod2)), N_Pod);
    matMult(1., A.row(idxPod1, 3*(countPod1+countPod2)).trans(), l.row(idxPod1, 3*(countPod1+countPod2)), n_Pod);
    lPl_Pod      += quadsum(l.row(idxPod1, 3*(countPod1+countPod2)));
    obsCountPod  += 3*(countPod1+countPod2);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
