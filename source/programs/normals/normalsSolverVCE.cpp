/***********************************************/
/**
* @file normalsSolverVCE.cpp
*
* @brief solve normal equations.
* Using cholesky decomposition.
* Relative weighting of different normal equations by variance component estimation (VCE).
*
* @author Torsten Mayer-Guerr
* @date 2003-03-26
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program accumulates \configClass{normalEquation}{normalEquationType}
and solves the total combined system.
The relative weigthing between the indivdual normals is determined iteratively
by means of variance component estimation (VCE). For a detailed description
of the used algorithm see \configClass{normalEquation}{normalEquationType}.

Besides the estimated parameter vector (\configFile{outputfileSolution}{matrix}) the
estimated accuracies (\configFile{outputfileSigmax}{matrix}) and the full covariance matrix
(\configFile{outputfileCovariance}{matrix}) can be saved. Also the combined normal system
can be written to \configFile{outputfileNormalEquation}{normalEquation}.

The \configFile{outputfileContribution}{matrix} is a matrix with rows for each estimated
parameter and columns for each \configClass{normalEquation}{normalEquationType}
and indicates the contribution of the individual normals to the estimated parameters.
Each row sum up to one.

See also \program{NormalsBuild}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileMatrix.h"
#include "classes/normalEquation/normalEquation.h"

/***** CLASS ***********************************/

/** @brief Solve normal equations.
* Using cholesky decomposition.
* relative weighting of different normal equations by variance component estimation (VCE).
* @ingroup programsGroup */
class NormalsSolverVCE
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(NormalsSolverVCE, PARALLEL, "solve normal equations (relative weighting by variance component estimation (VCE))", NormalEquation)

/***********************************************/

void NormalsSolverVCE::run(Config &config)
{
  try
  {
    FileName          fileNameSolution, fileNameSigmax;
    FileName          fileNameCovariance, fileNameNormals, fileNameContribution, fileNameVarianceFactors;
    NormalEquationPtr normals;
    FileName          fileNameX0;
    UInt              rhsNo;
    UInt              maxIter;
    UInt              blockSize;

    renameDeprecatedConfig(config, "outputfileNormalequation", "outputfileNormalEquation", date2time(2020, 6, 3));
    renameDeprecatedConfig(config, "normalequation",           "normalEquation",           date2time(2020, 6, 3));

    readConfig(config, "outputfileSolution",        fileNameSolution,        Config::OPTIONAL, "",     "parameter vector");
    readConfig(config, "outputfileSigmax",          fileNameSigmax,          Config::OPTIONAL, "",     "standard deviations of the parameters (sqrt of the diagonal of the inverse normal equation)");
    readConfig(config, "outputfileCovariance",      fileNameCovariance,      Config::OPTIONAL, "",     "full covariance matrix");
    readConfig(config, "outputfileContribution",    fileNameContribution,    Config::OPTIONAL, "",     "contribution of normal system components to the solution vector");
    readConfig(config, "outputfileVarianceFactors", fileNameVarianceFactors, Config::OPTIONAL, "",     "estimated variance factors as vector");
    readConfig(config, "outputfileNormalEquation",  fileNameNormals,         Config::OPTIONAL, "",     "the combined normal equation system");
    readConfig(config, "normalEquation",            normals,                 Config::MUSTSET,  "",     "");
    readConfig(config, "inputfileApproxSolution",   fileNameX0,              Config::OPTIONAL, "",     "to accelerate convergence");
    readConfig(config, "rightHandSideNumberVCE",    rhsNo,                   Config::DEFAULT,  "0",    "the right hand side number for estimation of variance factors");
    readConfig(config, "normalsBlockSize",          blockSize,               Config::DEFAULT,  "2048", "block size for distributing the normal equations, 0: one block");
    readConfig(config, "maxIterationCount",         maxIter,                 Config::DEFAULT,  "20",   "maximum number of iterations for variance component estimation");
    if(isCreateSchema(config)) return;

    logStatus<<"init normal equations"<<Log::endl;
    normals->init(blockSize);
    logInfo<<"  number of unknown parameters: "<<normals->parameterCount()<<Log::endl;
    logInfo<<"  number of right hand sides:   "<<normals->rightHandSideCount()<<Log::endl;

    // set approximate solution
    // ------------------------
    if(!fileNameX0.empty())
    {
      logStatus<<"read approximate solution <"<<fileNameX0<<">"<<Log::endl;
      Matrix x0;
      readFileMatrix(fileNameX0, x0);
      normals->setApproximateSolution(x0);
    }

    // Iteration
    // ---------
    Bool ready = FALSE;
    if(maxIter<1) maxIter = 1; // at least one step
    for(UInt iter=0; (iter<maxIter)&&(!ready); iter++)
    {
      logStatus<<iter+1<<". iteration step"<<Log::endl;

      logStatus<<"accumulate normal equations"<<Log::endl;
      ready = normals->build(rhsNo);
      logInfo<<"  observation count = "<<normals->observationCount()<<Log::endl;

      if(!fileNameNormals.empty())
      {
        logStatus<<"write normal equations to file <"<<fileNameNormals<<">"<<Log::endl;
        normals->write(fileNameNormals);
      }

      Vector sigmas = normals->varianceComponentFactors();
      if(Parallel::isMaster())
      {
        for (UInt i=0; i<sigmas.rows(); i++)
          sigmas(i) = std::sqrt(sigmas(i));
        constexpr UInt maxCount = 20;
        logInfo<<"  applied variance factors"<<Log::endl;
        for(UInt i=0; i<std::min(maxCount-2, sigmas.rows()); i++)
          logInfo<<"    "<<i+1<<": sigma = "<<sigmas(i)<<Log::endl;
        if(sigmas.rows()>maxCount)
          logInfo<<"    "<<maxCount-1<<": sigma = ..."<<Log::endl;
        if(sigmas.rows()>=maxCount-2)
        {
          logInfo<<"    "<<sigmas.rows()-1<<": sigma = "<<sigmas(sigmas.rows()-2)<<Log::endl;
          logInfo<<"    "<<sigmas.rows()  <<": sigma = "<<sigmas(sigmas.rows()-1)<<Log::endl;
        }

        if(!fileNameVarianceFactors.empty())
        {
          logStatus<<"write variance factors to <"<<fileNameVarianceFactors<<">"<<Log::endl;
          writeFileMatrix(fileNameVarianceFactors, sigmas);
        }
      }

      logStatus<<"solve normal equations"<<Log::endl;
      Matrix x = normals->solve();
      logInfo<<"  sigma (total) = "<<normals->aposterioriSigma()<<Log::endl;

      if(Parallel::isMaster() && !fileNameSolution.empty())
      {
        logStatus<<"write solution to <"<<fileNameSolution<<">"<<Log::endl;
        writeFileMatrix(fileNameSolution, x);
      }
    } // for(iter)

    // Inverse
    if(!fileNameSigmax.empty())
    {
      logStatus<<"inverte cholesky matrix and write standard deviations to <"<<fileNameSigmax<<">"<<Log::endl;
      Vector sigmax = normals->sigmaParameter();
      if(Parallel::isMaster())
        writeFileMatrix(fileNameSigmax, sigmax);
    }

    // Covariance matrix
    if(!fileNameCovariance.empty())
    {
      logStatus<<"compute covariance matrix and write to <"<<fileNameCovariance<<">"<<Log::endl;
      normals->writeCovariance(fileNameCovariance);
    }

    // contributions
    if(!fileNameContribution.empty())
    {
      logStatus<<"compute contributions and write to <"<<fileNameContribution<<">"<<Log::endl;
      Matrix contrib = normals->contribution();
      if(Parallel::isMaster())
        writeFileMatrix(fileNameContribution, contrib);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
