/***********************************************/
/**
* @file covarianceSst.cpp
*
* @brief Covariance matrix of satellite to satellite tracking observations.
*
* @author Torsten Mayer-Guerr
* @date 2010-07-18
*
*/
/***********************************************/

#define DOCSTRING_CovarianceSst

#include "base/import.h"
#include "config/configRegister.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"
#include "covarianceSst.h"

/***********************************************/

GROOPS_REGISTER_CLASS_WITHOUT_SUBS(CovarianceSst, "covarianceSstType")
GROOPS_READCONFIG_CLASS(CovarianceSst, "covarianceSstType")

/***********************************************/

CovarianceSst::CovarianceSst(Config &config, const std::string &name)
{
  try
  {
    FileName fileNameSigmaArc;
    FileName fileNameSigmaEpoch;
    FileName fileNameCovFunc;
    FileName fileNameCovMatrixSigmas;

    readConfigSequence(config, name, Config::MUSTSET, "", "");
    readConfig(config, "sigma",                        sigma,                     Config::DEFAULT,  "1", "general variance factor");
    readConfig(config, "inputfileSigmasPerArc",        fileNameSigmaArc,          Config::OPTIONAL, "",  "different accuaries for each arc (multplicated with sigma)");
    readConfig(config, "inputfileSigmasPerEpoch",      fileNameSigmaEpoch,        Config::OPTIONAL, "",  "different accuaries for each epoch (added)");
    readConfig(config, "inputfileCovarianceFunction",  fileNameCovFunc,           Config::OPTIONAL, "",  "covariance function in time");
    readConfig(config, "inputfileCovarianceMatrixArc", fileNamesCovarianceMatrix, Config::OPTIONAL, "",  "one matrix file per arc. Use {arcNo} as template");
    readConfig(config, "sigmasCovarianceMatrixArc",    fileNameCovMatrixSigmas,   Config::OPTIONAL, "",  "vector with one sigma for each covarianceMatrixArc");
    endSequence(config);
    if(isCreateSchema(config)) return;

    if(!fileNameSigmaArc.empty())
      readFileMatrix(fileNameSigmaArc, sigmaArc);

    if(!fileNameSigmaEpoch.empty())
      fileSigmaEpoch.open(fileNameSigmaEpoch);

    if(!fileNameCovFunc.empty())
      readFileMatrix(fileNameCovFunc, covFunction);

    if(!fileNameCovMatrixSigmas.empty())
    {
      readFileMatrix(fileNameCovMatrixSigmas, covMatrixSigmas);
      if(covMatrixSigmas.rows() != fileNamesCovarianceMatrix.size())
        throw(Exception("Number of sigmas not compatible with number of given arc-wise SST covariance matrices"));
    }
    else if(!fileNamesCovarianceMatrix.empty())
    {
      covMatrixSigmas = Vector(fileNamesCovarianceMatrix.size());
      for(UInt i=0; i<covMatrixSigmas.rows(); i++)
        covMatrixSigmas(i) = 1.0;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void CovarianceSst::testInput(const std::vector<Time> &times, const ObservationSigmaArc &sigmaEpoch, const_MatrixSliceRef covFunction, const_MatrixSliceRef W)
{
  if(sigmaEpoch.size() && (sigmaEpoch.size() != times.size()))
    throw(Exception("sigma per epoch not compatible with this arc number"));
  if(covFunction.size() && (covFunction.rows()<times.size()))
    throw(Exception("covariance function to short for this arc"));
  if(W.size() && ((W.rows() != times.size()) || (W.columns() != times.size())) )
    throw(Exception("covariance matrix W wrong size"));
}

/***********************************************/

Matrix CovarianceSst::covariance(UInt arcNo, const std::vector<Time> &times)
{
  try
  {
    Matrix C(times.size(), Matrix::SYMMETRIC, Matrix::UPPER);

    for(UInt i=0; i<fileNamesCovarianceMatrix.size(); i++)
    {
      Matrix M;
      FileName thisFile = fileNamesCovarianceMatrix.at(i).appendBaseName(".arc"+arcNo%"%03i"s);
      readFileMatrix(thisFile, M);
      if((M.rows() != C.rows()) || (M.columns() != C.columns()))
        throw(Exception("Arc-Wise SST covariance matrix <"+thisFile.str()+"> wrong size for arc."));
      axpy(std::pow(covMatrixSigmas(i), 2), M, C);
    }

    return covariance(times, (sigmaArc.size() ? sigmaArc(arcNo) : 1.), fileSigmaEpoch.readArc(arcNo), covFunction, C);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void CovarianceSst::decorrelate(UInt arcNo, const std::vector<Time> &times, const std::list<MatrixSlice> &A)
{
  try
  {
    Double weight = 1./this->sigma;
    if(sigmaArc.size() != 0)
      weight *= 1./this->sigmaArc(arcNo);

    // special case 1: diagonal matrix
    // -------------------------------
    if(fileNamesCovarianceMatrix.empty() && !covFunction.size() && !fileSigmaEpoch.arcCount())
    {
      for(MatrixSliceRef WA : A)
        if(WA.size()) WA *= weight;
      return;
    }

    // special case 2: diagonal matrix
    // -------------------------------
    if(fileNamesCovarianceMatrix.empty() && !covFunction.size())
    {
      ObservationSigmaArc sigmaEpoch = fileSigmaEpoch.readArc(arcNo);
      for(UInt i=0; i<sigmaEpoch.size(); i++)
        for(MatrixSliceRef WA : A)
          if(WA.size()) WA.row(i) *= weight/sigmaEpoch.at(i).sigma;
      return;
    }

    Matrix W = covariance(arcNo, times);
    cholesky(W);

    for(MatrixSliceRef WA : A)
      if(WA.size())
        triangularSolve(1., W.trans(), WA);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

MatrixSliceRef CovarianceSst::covariance(const std::vector<Time> &times, Double sigmaArc, const ObservationSigmaArc &sigmaEpoch,
                                         const_MatrixSliceRef covFunction, Matrix &Cov)
{
  try
  {
    testInput(times, sigmaEpoch, covFunction, Cov);

    if(!Cov.size())
      Cov = Matrix(times.size(), Matrix::SYMMETRIC, Matrix::UPPER);

    const Double sigma2 = pow(sigmaArc, 2);
    if(covFunction.size())
      for(UInt z=0; z<times.size(); z++)
        for(UInt s=z; s<times.size(); s++)
          Cov(z,s) += sigma2 * covFunction(s-z,1);

    for(UInt i=0; i<sigmaEpoch.size(); i++)
      Cov(i,i) += pow(sigmaEpoch.at(i).sigma,2);

    return Cov;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

MatrixSliceRef CovarianceSst::decorrelate(const std::vector<Time> &times, Double sigmaArc, const ObservationSigmaArc &sigmaEpoch,
                                          const_MatrixSliceRef covFunction, Matrix &W, const std::list<MatrixSlice> &A)
{
  try
  {
    covariance(times, sigmaArc, sigmaEpoch, covFunction, W);
    cholesky(W);

    for(MatrixSliceRef WA : A)
      if(WA.size())
        triangularSolve(1., W.trans(), WA);

    return W;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
