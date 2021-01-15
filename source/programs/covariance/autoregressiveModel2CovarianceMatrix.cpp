/***********************************************/
/**
* @file autoregressiveModel2CovarianceMatrix.cpp
*
* @brief Compute the covariance structure of a sequence of VAR(0) to VAR(p) models
*
* @author Andreas Kvas
* @date 2015-10-08
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program computes the covariance structure of a random process represented by an AR model sequence.
The covariance matrix is determined by accumulating the normal equations of all AR models in \config{autoregressiveModelSequence}
and inverting the combined normal equation matrix.
For each output file in \configFile{outputfileCovarianceMatrix}{matrix},
the covariance matrix of appropriate time lag is saved (the first file contains the auto-covariance,
second file cross covariance and so on). The matrix for lag $h$ describes the covariance between $x_{t-h}$ and $x_{t}$, i.e. $\Sigma(t-h, t)$.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileMatrix.h"
#include "parallel/matrixDistributed.h"
#include "misc/kalmanProcessing.h"

/***** CLASS ***********************************/

/** @brief Compute the covariance structure of a sequence of VAR(0) to VAR(p) models
* @ingroup programsGroup */
class AutoregressiveModel2CovarianceMatrix
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(AutoregressiveModel2CovarianceMatrix, SINGLEPROCESS, "Compute the covariance structure of a sequence of VAR(0) to VAR(p) models", Covariance)

/***********************************************/

void AutoregressiveModel2CovarianceMatrix::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    std::vector<FileName> fileNameOut;
    AutoregressiveModelSequencePtr arSequence;

    readConfig(config, "outputfileCovarianceMatrix",     fileNameOut,       Config::MUSTSET, "", "covariance matrix for each lag");
    readConfig(config, "autoregressiveModelSequence",    arSequence,        Config::MUSTSET, "", "AR model sequence");
    if(isCreateSchema(config)) return;

    logStatus<<"read autoregressive model sequence"<<Log::endl;
    const UInt blockCount = arSequence->maximumOrder()+1;
    const UInt dimension  = arSequence->dimension();
    auto blockIndex = MatrixDistributed::computeBlockIndex(dimension*blockCount, dimension);
    MatrixDistributed N(blockIndex, Parallel::selfCommunicator());

    logStatus<<"set up normal equations"<<Log::endl;
    for(UInt r = 0; r<N.blockCount(); r++)
      for(UInt c = r; c<N.blockCount(); c++)
        arSequence->distributedNormalsBlock(N.blockCount(), r, c, N.N(r, c));

    logStatus<<"compute inverse of normal equations"<<Log::endl;
    N.cholesky(FALSE/*timing*/);
    N.choleskyInverse(FALSE/*timing*/);
    N.choleskyProduct(FALSE/*timing*/);
    for(UInt column = 0; column<std::min(N.blockCount(), fileNameOut.size()); column++)
    {
      logStatus<<"write covariance matrix to <"<<fileNameOut.at(column)<<">"<<Log::endl;
      writeFileMatrix(fileNameOut.at(column), N.N(0, column));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
