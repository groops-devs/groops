/***********************************************/
/**
* @file noiseNormalsSolution.cpp
*
* @brief Create correlated errors from normal equations.
*
* @author Andreas Kvas
* @date 2017-06-27
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
The inverse of the normal matrix of \configFile{inputfileNormalEquation}{normalEquation}
represents the covariance matrix of the estimated parameters. This program generates a noise vector with
\begin{equation}
\M\Sigma(\M e) = \M N^{-1},
\end{equation}
if generated input noise is standard white noise.

The noise vector is computed with
\begin{equation}
\M e = \M W^{-T} \M z,
\end{equation}
where $\M z$ is the generated \configClass{noise}{noiseGeneratorType} and
$\M W$ is the cholesky upper triangle matrix of the normal matrix $\M N=\M W^T\M W$.
)";

/***********************************************/

#include "programs/program.h"
#include "parallel/matrixDistributed.h"
#include "files/fileMatrix.h"
#include "files/fileNormalEquation.h"
#include "classes/noiseGenerator/noiseGenerator.h"

/***** CLASS ***********************************/

/** * @brief Create a correlated error vector from a system of normal equations.
* @ingroup programsGroup */
class NoiseNormalsSolution
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(NoiseNormalsSolution, PARALLEL, "Create a correlated error vector from a system of normal equations.", Simulation, Noise, NormalEquation)
GROOPS_RENAMED_PROGRAM(NormalsSimulateNoise, NoiseNormalsSolution, date2time(2020, 7, 20))

/***********************************************/

void NoiseNormalsSolution::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName          outName, normalsName;
    NoiseGeneratorPtr noise;
    UInt              sampleCount;
    Bool              useEVD;

    renameDeprecatedConfig(config, "inputfileNormalequation",  "inputfileNormalEquation",  date2time(2020, 6, 3));

    readConfig(config, "outputfileNoise",         outName,     Config::MUSTSET,  "",  "generated noise as matrix: parameterCount x sampleCount");
    readConfig(config, "inputfileNormalEquation", normalsName, Config::MUSTSET,  "",  "");
    readConfig(config, "noise",                   noise,       Config::MUSTSET,  "",  "");
    readConfig(config, "sampleCount",             sampleCount, Config::DEFAULT,  "1", "number of samples to be generated");
    readConfig(config, "useEigenvalueDecomposition", useEVD,   Config::DEFAULT,  "0", "use eigenvalue decomposition");
    if(isCreateSchema(config)) return;

    if(useEVD)
    {
      if(!Parallel::isMaster(comm))
        return;

      Matrix N;
      Matrix n;
      NormalEquationInfo info;
      readFileNormalEquation(normalsName, info, N, n);
      logInfo<<"  number of parameters: "<<N.rows()<<Log::endl;

      Matrix rhs = noise->noise(n.rows(), sampleCount);

      Vector eig = eigenValueDecomposition(N);
      const Double maxE = maxabs(eig);

      Matrix y = N*rhs;
      for(UInt k = 0; k<y.rows(); k++)
      {
        Double factor = (eig(k)<(maxE*1e-13)) ? 0.0 : 1.0/std::sqrt(eig(k));
        y.row(k)*=factor;
      }
      Matrix out = N.trans()*y;

      logStatus<<"write noise vectors to file <"<<outName<<">"<<Log::endl;
      writeFileMatrix(outName, out);
      return;
    }

    // read normal equations
    // ---------------------
    logStatus<<"read normal equations <"<<normalsName<<">"<<Log::endl;
    MatrixDistributed normal;
    Matrix n;
    NormalEquationInfo info;
    readFileNormalEquation(normalsName, info, normal, n, comm);
    logInfo<<"  number of parameters: "<<normal.parameterCount()<<Log::endl;

    Matrix rhs = noise->noise(n.rows(), sampleCount);

    // compute cholesky
    // ----------------
    logStatus<<"compute cholesky decomposition"<<Log::endl;
    normal.cholesky(TRUE/*timing*/); // N = W^T W
    normal.triangularSolve(rhs);     // W^-1 e

    // save to file
    // ------------
    if(Parallel::isMaster(comm))
    {
      logStatus<<"write noise vectors to file <"<<outName<<">"<<Log::endl;
      writeFileMatrix(outName, rhs);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
