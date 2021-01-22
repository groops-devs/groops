/***********************************************/
/**
* @file covariancePod.cpp
*
* @brief Covariance matrix of kinematic orbits.
*
* @author Torsten Mayer-Guerr
* @date 2010-07-18
*
*/
/***********************************************/

#define DOCSTRING_CovariancePod

#include "base/import.h"
#include "config/configRegister.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"
#include "covariancePod.h"

/***********************************************/

GROOPS_REGISTER_CLASS_WITHOUT_SUBS(CovariancePod, "covariancePodType")
GROOPS_READCONFIG_CLASS(CovariancePod, "covariancePodType")

/***********************************************/

CovariancePod::CovariancePod(Config &config, const std::string &name)
{
  try
  {
    FileName fileNameSigmaArc, fileNameSigmaEpoch, fileNameCovPodEpoch, fileNameCovFunc;

    readConfigSequence(config, name, Config::MUSTSET, "", "");
    readConfig(config, "sigma",                        sigma,               Config::DEFAULT,  "1", "general variance factor");
    readConfig(config, "inputfileSigmasPerArc",        fileNameSigmaArc,    Config::OPTIONAL, "",  "different accuaries for each arc (multplicated with sigma)");
    readConfig(config, "inputfileSigmasPerEpoch",      fileNameSigmaEpoch,  Config::OPTIONAL, "",  "different accuaries for each epoch (added)");
    readConfig(config, "inputfileCovarianceFunction",  fileNameCovFunc,     Config::OPTIONAL, "",  "covariances in time for along, cross, and radial direction");
    readConfig(config, "inputfileCovariancePodEpoch",  fileNameCovPodEpoch, Config::OPTIONAL, "",  "3x3 epoch wise covariances");
    endSequence(config);
    if(isCreateSchema(config)) return;

    if(!fileNameSigmaArc.empty())
      readFileMatrix(fileNameSigmaArc, sigmaArc);

    if(!fileNameSigmaEpoch.empty())
      fileSigmaEpoch.open(fileNameSigmaEpoch);

    if(!fileNameCovPodEpoch.empty())
      fileCovPodEpoch.open(fileNameCovPodEpoch);

    if(!fileNameCovFunc.empty())
      readFileMatrix(fileNameCovFunc, covFunction);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void CovariancePod::testInput(const OrbitArc &pod, const ObservationSigmaArc &sigmaEpoch, const Covariance3dArc &covPod, const_MatrixSliceRef covFunction)
{
  if(sigmaEpoch.size() && (sigmaEpoch.size() != pod.size()))
    throw(Exception("sigma per epoch not compatible with this arc number"));
  if((covPod.size() != 0) && (covPod.size() != pod.size()))
    throw(Exception("orbit and CovariancePodEpoch are not compatible"));
  if(covFunction.size() && (covFunction.rows()<pod.size()))
    throw(Exception("covariance function to short for this arc"));
}

/***********************************************/

Matrix CovariancePod::covariance(UInt arcNo, const OrbitArc &pod)
{
  try
  {
    ObservationSigmaArc sigmaEpoch = fileSigmaEpoch.readArc(arcNo);
    Covariance3dArc     covPod     = fileCovPodEpoch.readArc(arcNo);
    testInput(pod, sigmaEpoch, covPod, covFunction);
    if(sigmaArc.size() && (arcNo>=sigmaArc.size()))
      throw(Exception("sigmasPerArc contain not enough rows for this arc number"));

    Double sigma2 = pow(this->sigma, 2);
    if(sigmaArc.size())
      sigma2 *= pow(this->sigmaArc(arcNo), 2);

    Matrix C(3*pod.size(), Matrix::TRIANGULAR);
    for(UInt i=0; i<sigmaEpoch.size(); i++)
      for(UInt k=0; k<3; k++)
        C(3*i+k,3*i+k) += pow(sigmaEpoch.at(i).sigma, 2);

    // special case 1: diagonal matrix
    // -------------------------------
    if((covPod.size()==0) && (covFunction.size()==0))
    {
      for(UInt i=0; i<C.rows(); i++)
        C(i,i) += sigma2;
      return C;
    }

    // special case 2: epoch block diagonal matrix
    // -------------------------------------------
    if(covPod.size() && (covFunction.size()==0))
    {
      for(UInt i=0; i<pod.size(); i++)
      {
        Matrix C3x3 = covPod.at(i).covariance.matrix();
        axpy(sigma2, C3x3, C.slice(3*i,3*i,3,3));
      }
      return C;
    }

    // covariance function in orbit system
    // -----------------------------------
    const Double sampling = covFunction(1,0)-covFunction(0,0);
    for(UInt z=0; z<pod.size(); z++)
      for(UInt s=z; s<pod.size(); s++)
      {
        UInt idx = static_cast<UInt>(round((pod.at(s).time-pod.at(z).time).seconds()/sampling));
        C(3*z+0, 3*s+0) += sigma2 * covFunction(idx, 1+0);
        C(3*z+1, 3*s+1) += sigma2 * covFunction(idx, 1+1);
        C(3*z+2, 3*s+2) += sigma2 * covFunction(idx, 1+2);
      }
    fillSymmetric(C);

    // rotate and decorrelate epoch wise residuals
    // -------------------------------------------
    for(UInt i=0; i<pod.size(); i++)
    {
      // orbit system: z: radial, x: along, y: cross
      Vector3d x;
      if(i==0)
        x = pod.at(i+1).position - pod.at(i).position;
      else
        x = pod.at(i).position - pod.at(i-1).position;
      Vector3d z = normalize(pod.at(i).position);
      Vector3d y = normalize(crossProduct(z, x));
      x = crossProduct(y, z);
      Matrix D = Rotary3d(x,y).matrix().trans(); // Rot from CRF to SRF

      // epoch wise decorrelation
      if(covPod.size()!=0)
      {
        Matrix V = covPod.at(i).covariance.matrix();    // 3x3 Epoch covariance
        V.setType(Matrix::SYMMETRIC, Matrix::UPPER);
        Vector eigen = eigenValueDecomposition(V);
        Matrix V1 = V;
        V1.column(0) *= std::sqrt(eigen(0));
        V1.column(1) *= std::sqrt(eigen(1));
        V1.column(2) *= std::sqrt(eigen(2));
        D = D*V1*V.trans();           // W^-T im LORF
      }

      copy(D.trans() * C.row(3*i,3), C.row(3*i,3));
      copy(C.column(3*i,3) * D,      C.column(3*i,3));
    }

    return C;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void CovariancePod::decorrelate(UInt arcNo, const OrbitArc &pod, const std::list<MatrixSlice> &A)
{
  try
  {
    ObservationSigmaArc sigmaEpoch = fileSigmaEpoch.readArc(arcNo);
    Covariance3dArc     covPod     = fileCovPodEpoch.readArc(arcNo);
    testInput(pod, sigmaEpoch, covPod, covFunction);
    if(sigmaArc.size() && (arcNo>=sigmaArc.size()))
      throw(Exception("sigmasPerArc contain not enough rows for this arc number"));

    Double weight = 1./this->sigma;
    if(sigmaArc.size() != 0)
      weight *= 1./this->sigmaArc(arcNo);

    // special case 1: diagonal matrix
    // -------------------------------
    if((covPod.size()==0) && (sigmaEpoch.size()==0) && (covFunction.size()==0))
    {
      for(MatrixSliceRef WA : A)
        if(WA.size()) WA *= weight;
      return;
    }

    // special case 2: epoch block diagonal matrix
    // -------------------------------------------
    if(covPod.size() && (covFunction.size()==0))
    {
      for(UInt i=0; i<pod.size(); i++)
      {
        Matrix W = covPod.at(i).covariance.matrix();
        W.setType(Matrix::SYMMETRIC);

        if(sigmaEpoch.size())
        {
          W(0,0) += sigmaEpoch.at(i).sigma;
          W(1,1) += sigmaEpoch.at(i).sigma;
          W(2,2) += sigmaEpoch.at(i).sigma;
        }

        cholesky(W);
        for(MatrixSliceRef WA : A)
          if(WA.size())
            triangularSolve(weight, W.trans(), WA.row(3*i,3));
      }
      return;
    }

    decorrelate(pod, 1/weight, sigmaEpoch, covPod, covFunction, A);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix CovariancePod::decorrelate(const OrbitArc &pod, Double sigmaArc, const ObservationSigmaArc &sigmaEpoch,
                                  const Covariance3dArc &covPod, const_MatrixSliceRef covFunction,
                                  const std::list<MatrixSlice> &A)
{
  try
  {
    testInput(pod, sigmaEpoch, covPod, covFunction);

    // rotate and decorrelate epoch wise residuals
    // -------------------------------------------
    for(UInt i=0; i<pod.size(); i++)
    {
      // orbit system: z: radial, x: along, y: cross
      Vector3d x;
      if(i==0)
        x = pod.at(i+1).position - pod.at(i).position;
      else
        x = pod.at(i).position - pod.at(i-1).position;
      Vector3d z = normalize(pod.at(i).position);
      Vector3d y = normalize(crossProduct(z, x));
      x = crossProduct(y, z);
      Matrix D = Rotary3d(x,y).matrix().trans(); // Rot from CRF to SRF

      // epoch wise decorrelation
      if(covPod.size()!=0)
      {
        Matrix V = covPod.at(i).covariance.matrix(); // 3x3 Epoch covariance
        V.setType(Matrix::SYMMETRIC, Matrix::UPPER);
        Vector eigen = eigenValueDecomposition(V);
        Matrix V1 = V;
        V1.column(0) *= 1./sqrt(eigen(0));
        V1.column(1) *= 1./sqrt(eigen(1));
        V1.column(2) *= 1./sqrt(eigen(2));
        D = D*V1*V.trans();           // W^-T im LORF
      }

      for(MatrixSliceRef WA : A)
        if(WA.size())
          copy(D * WA.row(3*i,3), WA.row(3*i,3));
    }

    // covariance function in orbit system
    // -----------------------------------
    Matrix W(3*pod.size(), Matrix::SYMMETRIC, Matrix::UPPER);
    if(covFunction.size())
    {
      const Double sampling = covFunction(1,0)-covFunction(0,0);
      for(UInt z=0; z<pod.size(); z++)
        for(UInt s=z; s<pod.size(); s++)
        {
          UInt idx = static_cast<UInt>(round((pod.at(s).time-pod.at(z).time).seconds()/sampling));
          W(3*z+0, 3*s+0) = sigmaArc*sigmaArc * covFunction(idx, 1+0);
          W(3*z+1, 3*s+1) = sigmaArc*sigmaArc * covFunction(idx, 1+1);
          W(3*z+2, 3*s+2) = sigmaArc*sigmaArc * covFunction(idx, 1+2);
        }
    }

    for(UInt i=0; i<sigmaEpoch.size(); i++)
      for(UInt k=0; k<3; k++)
        W(3*i+k,3*i+k) += pow(sigmaEpoch.at(i).sigma, 2);

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
