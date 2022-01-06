/***********************************************/
/**
* @file preprocessingVariationalEquationOrbitFit.cpp
*
* @brief Fit variational equations to orbit observations.
*
* @author Torsten Mayer-Guerr
* @date 2012-05-30
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program fits an \configFile{inputfileVariational}{variationalEquation} to an observed \configFile{inputfileOrbit}{instrument} by estimating parameters
in a least squares adjustment. Additional to the initial satellite state for each arc, these parameters can be
\configClass{parametrizationGravity}{parametrizationGravityType}, satellite \configClass{parametrizationAcceleration}{parametrizationAccelerationType}
and stochastic pulses (velocity jumps) at given times, \configClass{stochasticPulse}{timeSeriesType}. The estimated parameters can be stored with
\configFile{outputfileSolution}{matrix} and an extra file with the parameter names is created. The fitted orbit is written
as new reference in \configFile{outputfileVariational}{variationalEquation} and additionally in \configFile{outputfileOrbit}{instrument}.

The observed orbit positions (\configFile{inputfileOrbit}{instrument}) together with the epoch wise covariance matrix
(\configFile{inputfileCovariancePodEpoch}{instrument}) must be splitted in the same arcs as the variational equations but not
necessarily uniform distributed (use irregularData in \program{InstrumentSynchronize}). An iterative downweighting of
outliers is performed by M-Huber method.

The observation equations (parameter sensitity matrix) are computed by integration of the variational equations
(\configFile{inputfileVariational}{variationalEquation}) using a polynomial with \config{integrationDegree} and interpolated to the
observation epochs using a polynomial with \config{interpolationDegree}.

All parameters used here must be reestimated in the full least squares adjustment
for the gravity field determination to get a solution which is not biased towards the reference field.
The solution of additional estimations are relative (deltas) as the parameters are already used as Taylor point
in the reference orbit.

See also \program{PreprocessingVariationalEquation}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/polynomial.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"
#include "files/fileVariationalEquation.h"
#include "files/fileParameterName.h"
#include "classes/ephemerides/ephemerides.h"
#include "classes/parametrizationGravity/parametrizationGravity.h"
#include "classes/parametrizationAcceleration/parametrizationAcceleration.h"
#include "classes/timeSeries/timeSeries.h"
#include "misc/varianceComponentEstimation.h"
#include "misc/observation/variationalEquationFromFile.h"

/***** CLASS ***********************************/

/** @brief Fit variational equations to orbit observations
* @ingroup programsGroup */
class PreprocessingVariationalEquationOrbitFit
{
public:
  VariationalEquationFromFile    variationalEquationFromFile;
  InstrumentFile                 podFile;
  InstrumentFile                 covPodEpochFile;
  EphemeridesPtr                 ephemerides;
  ParametrizationAccelerationPtr parameterAcceleration;
  ParametrizationGravityPtr      parameterGravity;
  UInt                           interpolationDegree;
  UInt                           arcCount;

  // normal equations
  // ----------------
  Matrix N;           // =A'PA, Normal matrix
  Vector n;           // =A'Pl, right hand side
  Double lPl;         // =l'Pl, weighted norm of the observations
  UInt   obsCount;    // number of observations
  UInt   outlierCount;
  Vector x;
  Double sigma0;

  void run(Config &config, Parallel::CommunicatorPtr comm);

  void buildNormals(UInt arcNo);
};

GROOPS_REGISTER_PROGRAM(PreprocessingVariationalEquationOrbitFit, PARALLEL, "fit variational equations to orbit observations", Preprocessing, VariationalEquation)

/***********************************************/
/***********************************************/

void PreprocessingVariationalEquationOrbitFit::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName fileNameOutVariational, fileNameOutOrbit, fileNameOutSolution;
    FileName fileNameInVariational;
    FileName podName, covPodEpochName;
    UInt              integrationDegree;
    UInt              iterCount;
    std::vector<Time> stochasticPulse;
    TimeSeriesPtr     stochasticPulsePtr;

    renameDeprecatedConfig(config, "representation", "parametrizationGravity",      date2time(2020, 6, 3));
    renameDeprecatedConfig(config, "parameter",      "parametrizationAcceleration", date2time(2020, 6, 3));

    readConfig(config, "outputfileVariational",       fileNameOutVariational, Config::MUSTSET,  "",    "approximate position and integrated state matrix");
    readConfig(config, "outputfileOrbit",             fileNameOutOrbit,       Config::OPTIONAL, "",    "integrated orbit");
    readConfig(config, "outputfileSolution",          fileNameOutSolution,    Config::OPTIONAL, "",    "estimated calibration and state parameters");
    readConfig(config, "inputfileVariational",        fileNameInVariational,  Config::MUSTSET,  "",    "approximate position and integrated state matrix");
    readConfig(config, "inputfileOrbit",              podName,                Config::MUSTSET,  "",    "kinematic positions of satellite as observations");
    readConfig(config, "inputfileCovariancePodEpoch", covPodEpochName,        Config::OPTIONAL, "",    "3x3 epoch wise covariances");
    readConfig(config, "ephemerides",                 ephemerides,            Config::OPTIONAL, "jpl", "may be needed by parametrizationAcceleration");
    readConfig(config, "parametrizationGravity",      parameterGravity,       Config::DEFAULT,  "",    "gravity field parametrization");
    readConfig(config, "parametrizationAcceleration", parameterAcceleration,  Config::DEFAULT,  "",    "orbit force parameters");
    readConfig(config, "stochasticPulse",             stochasticPulsePtr,     Config::DEFAULT,  "",    "");
    readConfig(config, "integrationDegree",           integrationDegree,      Config::DEFAULT,  "7",   "integration of forces by polynomial approximation of degree n");
    readConfig(config, "interpolationDegree",         interpolationDegree,    Config::DEFAULT,  "7",   "orbit interpolation by polynomial approximation of degree n");
    readConfig(config, "iterationCount",              iterCount,              Config::DEFAULT,  "10",  "for the estimation of calibration parameter and error PSD");
    if(isCreateSchema(config)) return;

    if(integrationDegree%2 == 0)
      throw(Exception("polnomial degree for integration must be odd."));

    if(stochasticPulsePtr)
      stochasticPulse = stochasticPulsePtr->times();

    // init
    // ----
    podFile.open(podName);
    covPodEpochFile.open(covPodEpochName);
    InstrumentFile::checkArcCount({podFile, covPodEpochFile});
    variationalEquationFromFile.open(fileNameInVariational, parameterGravity, parameterAcceleration, stochasticPulse, ephemerides, integrationDegree);

    // =============================================

    x = Vector(variationalEquationFromFile.parameterCount());
    sigma0 = 1;
    for(UInt iter=0; iter<iterCount; iter++)
    {
      // build normals
      // -------------
      logStatus<<"accumulate normal equations"<<Log::endl;
      N            = Matrix(variationalEquationFromFile.parameterCount(), Matrix::SYMMETRIC);
      n            = Vector(variationalEquationFromFile.parameterCount());
      lPl          = 0;
      obsCount     = 0;
      outlierCount = 0;

      Parallel::forEach(podFile.arcCount(), [this](UInt arcNo) {buildNormals(arcNo);}, comm);

      Parallel::reduceSum(N,            0, comm);
      Parallel::reduceSum(n,            0, comm);
      Parallel::reduceSum(lPl,          0, comm);
      Parallel::reduceSum(obsCount,     0, comm);
      Parallel::reduceSum(outlierCount, 0, comm);

      // Estimate parameters
      // -------------------
      if(Parallel::isMaster(comm))
      {
        // Regularize unused parameters
        for(UInt i=0; i<N.rows(); i++)
          if(N(i,i) == 0)
            N(i,i) = 1.0;

        logStatus<<"solve system of normal equations"<<Log::endl;
        x = solve(N,n);
        sigma0 = Vce::standardDeviation(lPl-inner(n,x), obsCount-x.rows(), 2.5/*huber*/, 1./*huberPower*/);
        logInfo<<"  aposteriori sigma: "<<sigma0<<Log::endl;
        logInfo<<"  outlier "<<outlierCount<<" of "<<obsCount<<" ("<<100.*outlierCount/obsCount<<"%)"<<Log::endl;
      }

      Parallel::broadCast(outlierCount, 0, comm);
      if((iter>0) && (outlierCount==0))
        break;

      logInfo<<"  parameter count = "<<variationalEquationFromFile.parameterCount()<<Log::endl;
      Parallel::broadCast(x, 0, comm);
      Parallel::broadCast(sigma0, 0, comm);
    } // for(iter)

    if(Parallel::isMaster(comm) && !fileNameOutSolution.empty())
    {
      logStatus<<"write solution to <"<<fileNameOutSolution<<">"<<Log::endl;
      writeFileMatrix(fileNameOutSolution, x);

      std::vector<ParameterName> parameterName;
      variationalEquationFromFile.parameterName(parameterName);
      writeFileParameterName(fileNameOutSolution.replaceFullExtension(".parameterName.txt"), parameterName);
    }

    // =============================================================================

    // Improve orbits with estimated parameters
    // ----------------------------------------
    std::vector<VariationalEquationArc> arcs(variationalEquationFromFile.arcCount());
    Parallel::forEach(arcs, [this](UInt arcNo) {return variationalEquationFromFile.refineVariationalEquationArc(arcNo, x);}, comm);

    // =============================================================================

    if(Parallel::isMaster(comm) && !fileNameOutVariational.empty())
    {
      logStatus<<"write variational equation to file <"<<fileNameOutVariational<<">"<<Log::endl;
      writeFileVariationalEquation(fileNameOutVariational, variationalEquationFromFile.satellite(), arcs);
    }

    // =============================================

    if(Parallel::isMaster(comm) && !fileNameOutOrbit.empty())
    {
      logStatus<<"write orbit to file <"<<fileNameOutOrbit<<">"<<Log::endl;
      std::list<Arc> arcList;
      for(UInt arcNo=0; arcNo<arcs.size(); arcNo++)
        arcList.push_back( arcs.at(arcNo).orbitArc() );
      InstrumentFile::write(fileNameOutOrbit, arcList);
    }

    // =============================================
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void PreprocessingVariationalEquationOrbitFit::buildNormals(UInt arcNo)
{
  try
  {
    OrbitArc pod = podFile.readArc(arcNo);
    if(pod.size() == 0)
      return;

    Vector l(3*pod.size());
    for(UInt k=0; k<pod.size(); k++)
    {
      l(3*k+0) = pod.at(k).position.x();
      l(3*k+1) = pod.at(k).position.y();
      l(3*k+2) = pod.at(k).position.z();
    }

    std::vector<Time> timePod = pod.times();
    VariationalEquationFromFile::ObservationEquation eqn = variationalEquationFromFile.integrateArc(timePod.front(), timePod.back(), TRUE/*position*/, FALSE/*velocity*/);
    Polynomial polynomial(eqn.times, interpolationDegree);
    l -= polynomial.interpolate(timePod, eqn.pos0, 3); // reference orbit
    Matrix A = polynomial.interpolate(timePod, eqn.PosDesign, 3);

    // decorrelation
    Covariance3dArc covPod = covPodEpochFile.readArc(arcNo);
    Arc::checkSynchronized({pod, covPod});
    for(UInt i=0; i<covPod.size(); i++)
    {
      Matrix W = covPod.at(i).covariance.matrix();
      W.setType(Matrix::SYMMETRIC);
      cholesky(W);

      triangularSolve(1., W.trans(), l.row(3*i,3));
      triangularSolve(1., W.trans(), A.row(3*i,3));
    }

    // downweight outliers
    if(quadsum(x))
    {
      const Double huber = 2.5;
      Vector e = l;
      matMult(-1, A, x, e);
      for(UInt k=0; k<pod.size(); k++)
      {
        Double s = sqrt(quadsum(e.row(3*k,3))/3);
        if(s>huber*sigma0)
        {
          l.row(3*k,3) *= huber*sigma0/s;
          A.row(3*k,3) *= huber*sigma0/s;
          outlierCount += 3;
        }
      }
    }

    lPl += quadsum(l);
    obsCount += l.rows();
    rankKUpdate(1., A, N);
    matMult(1., A.trans(), l, n);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
