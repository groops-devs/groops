/***********************************************/
/**
* @file covarianceMatrix2Correlation.cpp
*
* @brief Compute the correlation matrix from a covariance matrix
*
* @author Andreas Kvas
* @date 2018-01-10
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program computes the pearson correlation coefficient
\begin{equation}
  \rho_{ij} = \frac{\sigma_{ij}}{\sigma_i \sigma_j}
\end{equation}
from a given covariance matrix stored in \configFile{inputfileCovarianceMatrix}{matrix}.
The result is stored in \configFile{outputfileCorrelationMatrix}{matrix}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileMatrix.h"

/***** CLASS ***********************************/

/** @brief Compute the correlation matrix from a covariance matrix.
* @ingroup programsGroup */
class CovarianceMatrix2Correlation
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(CovarianceMatrix2Correlation, PARALLEL, "compute the correlation matrix from a covariance matrix", Covariance)

/***********************************************/

void CovarianceMatrix2Correlation::run(Config &config)
{
  try
  {
    FileName fileNameOut, fileNameIn;

    readConfig(config, "outputfileCorrelationMatrix", fileNameOut, Config::MUSTSET, "",  "correlation matrix");
    readConfig(config, "inputfileCovarianceMatrix",   fileNameIn,  Config::MUSTSET, "",  "covariance matrix");
    if(isCreateSchema(config)) return;

    Matrix covarianceMatrix;
    readFileMatrix(fileNameIn, covarianceMatrix);
    if(covarianceMatrix.getType() != Matrix::SYMMETRIC)
      throw(Exception("Covariance matrix must be symmetric."));
    fillSymmetric(covarianceMatrix);

    Matrix correlationMatrix(covarianceMatrix.rows(), Matrix::SYMMETRIC, Matrix::UPPER);
    for(UInt r=0; r<covarianceMatrix.rows(); r++)
      for(UInt c=r; c<covarianceMatrix.columns(); c++)
        correlationMatrix(r,c) = covarianceMatrix(r,c)/std::sqrt(covarianceMatrix(r,r)*covarianceMatrix(c,c));

    logStatus<<"write correlation matrix to <"<<fileNameOut<<">"<<Log::endl;
    writeFileMatrix(fileNameOut, correlationMatrix);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
