/***********************************************/
/**
* @file CovarianceMatrix2AutoregressiveModel.cpp
*
* @brief Compute a VAR(p) model from covariance matrices
*
* @author Andreas Kvas
* @date 2015-10-19
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program computes a VAR(p) model from empirical covariance matrices.
The \configFile{inputfileCovarianceMatrix}{matrix} represent the covariance structure of the process:
the first file should contain the auto-covariance, the second the cross-covariance of lag one,
the next cross-covariance of lag two and so on.

Cross-covariance matrices $\Sigma_{\Delta_k}$ are defined as the cross-covariance between epoch $t-k$ and $t$.
If the process realizations $x_{t}$ are arrange by ascending time stamps
($\{\dots, x_{t-2}, x_{t-1}, x_{t}, x_{t+1}, x_{t+2},\dots\}$),
the covariance structure of the (stationary) process is therefore given by
\begin{equation}
\begin{bmatrix}
\Sigma & \Sigma_{\Delta_1} & \Sigma_{\Delta_2} & \cdots \\
\Sigma_{\Delta_1}^T & \Sigma & \Sigma_{\Delta_1} &  \cdots \\
\Sigma_{\Delta_2}^T & \Sigma_{\Delta_1}^T & \Sigma & \cdots \\
\vdots & \vdots & \vdots & \ddots \\
\end{bmatrix}.
\end{equation}

The estimate AR model is saved as single matrix \config{outputfileAutoregressiveModel} according to the GROOPS AR model conventions.
)";

/***********************************************/

#include "programs/program.h"
#include "parallel/matrixDistributed.h"
#include "files/fileMatrix.h"

/***** CLASS ***********************************/

/** @brief Compute least squares prediction matrices (B, Q) from covariance matrices
*
* @ingroup programsGroup */
class CovarianceMatrix2AutoregressiveModel
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(CovarianceMatrix2AutoregressiveModel, SINGLEPROCESS, "Compute a VAR(p) model from covariance matrices", Covariance)

/***********************************************/

void CovarianceMatrix2AutoregressiveModel::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    std::vector<FileName> inNameCovariance;
    FileName outNameModel;

    readConfig(config, "outputfileAutoregressiveModel",   outNameModel,      Config::MUSTSET,   "",  "coefficients and white noise covariance of AR(p) model");
    readConfig(config, "inputfileCovarianceMatrix",       inNameCovariance,  Config::MUSTSET,   "",  "file name of covariance matrix");
    if(isCreateSchema(config)) return;

    logStatus<<"read covariance matrices"<<Log::endl;
    const UInt order = inNameCovariance.size()-1;
    std::vector<Matrix> C(order+1);
    for(UInt k = 0; k<order+1; k++)
      readFileMatrix(inNameCovariance.at(k), C.at(k));

    if(C.front().getType() != Matrix::SYMMETRIC)
      throw(Exception("Auto-covariance matrix must be symmetric."));

    logStatus<<"compute autoregressive model of order <"<<order<<">"<<Log::endl;
    if(order == 0)
    {
      Matrix model = matrixSquareRootInverse(C.front());
      logStatus<<"write autoregressive model to <"<<outNameModel<<">"<<Log::endl;
      writeFileMatrix(outNameModel, model);
      return;
    }
    const UInt dim = C.front().rows();
    Matrix n(dim*order, dim); // right hand side
    MatrixDistributed N;
    auto blockIndex = MatrixDistributed::computeBlockIndex(dim*order, dim);
    N.initEmpty(blockIndex, Parallel::selfCommunicator());
    for(UInt r = 0; r<N.blockCount(); r++)
    {
      copy(C.at(r+1), n.row(r*dim, dim));
      for(UInt c = r; c<N.blockCount(); c++)
      {
        N.setBlock(r, c);
        N.N(r, c) = C.at(c-r).trans(); // block toeplitz
      }
    }
    Matrix X = N.solve(n, TRUE/*timing*/);

    logStatus<<"compute white noise covariance"<<Log::endl;
    Matrix Q = C.front();
    fillSymmetric(Q);
    Q.setType(Matrix::GENERAL);
    matMult(-1.0, X.trans(), n, Q);
    Q.setType(Matrix::SYMMETRIC);

    try
    {
      Matrix W = Q; cholesky(W);
    }
    catch(std::exception &/*e*/)
    {
      logWarning<<"white noise covariance is not positive definite"<<Log::endl;
    }

    Matrix F = matrixSquareRootInverse(Q); // F = Q^(-1/2)
    Matrix model(dim, (order+1)*dim); // [Ap, ..., A1, A0]
    copy(F, model.column(order*dim, dim));
    for(UInt k = 0; k<order; k++)
      matMult(-1.0, F, X.row((order-k-1)*dim, dim).trans(), model.column(k*dim, dim)); // Ak = -Q^(1/2)*B_k, B0 = -I

    logStatus<<"write autoregressive model to <"<<outNameModel<<">"<<Log::endl;
    writeFileMatrix(outNameModel, model);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

