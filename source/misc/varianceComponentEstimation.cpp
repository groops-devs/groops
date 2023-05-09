/***********************************************/
/**
* @file varianceComponentEstimation.cpp
*
* @brief Variance Component Estimation (VCE).
*
* @author Torsten Mayer-Guerr
* @date 2011-10-17
*
*/
/***********************************************/

#include <random>
#include "base/import.h"
#include "inputOutput/logging.h"
#include "files/fileMatrix.h"
#include "varianceComponentEstimation.h"

/***********************************************/

Matrix Vce::monteCarlo(UInt rows, UInt columns)
{
  // init random generator
  static std::mt19937_64 generator;
  static Bool init = FALSE;
  if(!init)
  {
   std::random_device randomDevice;
   generator.seed(randomDevice());
   init = TRUE;
  }

  Matrix e(rows, columns);
  std::uniform_int_distribution<Int> binary(0,1);
  for(UInt i=0; i<rows; i++)
    for(UInt k=0; k<columns; k++)
      e(i, k) = binary(generator) * 2 - 1;
  return 1./std::sqrt(columns) * e;
}

/***********************************************/

Double Vce::standardDeviation(Double ePe, Double redundancy, Double huber, Double huberPower)
{
  static Double huber_      = NAN;
  static Double huberPower_ = NAN;
  static Double sum         = NAN;
  if((huber_ == huber) && (huberPower_ == huberPower))
    return std::sqrt(ePe/redundancy/sum);

  constexpr Double dx = 1e-4;
  huber_      = huber;
  huberPower_ = huberPower;
  sum         = 0;
  Double x    = dx/2;

  // variance of normal distribution
  for(; x<std::min(huber, 10.); x+=dx)
    sum += x*x * std::exp(-0.5*x*x) * dx;
  // variance of downweighted normal distribution
  for(; x<10.; x+=dx)
    sum += x*x*std::pow(x/huber, -2*huberPower) * std::exp(-0.5*x*x) * dx;
  sum *= 2./std::sqrt(2*PI);
  return std::sqrt(ePe/redundancy/sum);
}

/***********************************************/

Matrix Vce::robustLeastSquares(const_MatrixSliceRef A, const_MatrixSliceRef l, UInt countGroup,
                               Double huber, Double huberPower, UInt maxIter, Vector &sigma)
{
  try
  {
    std::vector<UInt> indexGroup(l.rows()/countGroup+1);
    indexGroup.at(0) = 0;
    for(UInt i=1; i<indexGroup.size(); i++)
      indexGroup.at(i) = indexGroup.at(i-1) + countGroup;
    return robustLeastSquares(A, l, indexGroup, huber, huberPower, maxIter, sigma);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix Vce::robustLeastSquares(const_MatrixSliceRef A, const_MatrixSliceRef l, const std::vector<UInt> &indexGroup,
                               Double huber, Double huberPower, UInt maxIter, Vector &sigma)
{
  try
  {
    Matrix x;
    UInt   countOutlier = 0;
    Double sigma0 = 0;
    sigma = Vector(indexGroup.size()-1, 1.);
    for(UInt iter=0; iter<maxIter; iter++)
    {
      // weighting
      Matrix Wl = l;
      Matrix WA = A;
      for(UInt i=0; i<sigma.rows(); i++)
      {
        Wl.row(indexGroup.at(i), indexGroup.at(i+1)-indexGroup.at(i)) *= 1./sigma(i);
        WA.row(indexGroup.at(i), indexGroup.at(i+1)-indexGroup.at(i)) *= 1./sigma(i);
      }

      // QR decomposition
      const Vector tau = QR_decomposition(WA);
      QTransMult(WA, tau, Wl);           // transform observations: l:= Q'l
      x = Wl.row(0, WA.columns());
      triangularSolve(1., WA.row(0, WA.columns()), x);
      Wl.row(0, WA.columns()).setNull(); // residuals: remove WB*x
      QMult(WA, tau, Wl);                // back transformation
      generateQ(WA, tau);                // for redundancies

      if(sigma0 == 0.)
        sigma0 = std::sqrt(quadsum(Wl)/(Wl.size()-x.size()));

      // outlier detection
      UInt   countOutlierNew = 0;
      Double ePeSum = 0.;
      Double rSum   = 0.;
      for(UInt i=0; i<sigma.rows(); i++)
      {
        const Double ePe = quadsum(Wl.row(indexGroup.at(i), indexGroup.at(i+1)-indexGroup.at(i)))/Wl.columns();
        const Double r   = indexGroup.at(i+1)-indexGroup.at(i) - quadsum(WA.row(indexGroup.at(i), indexGroup.at(i+1)-indexGroup.at(i)));
        const Double s   = std::sqrt(ePe/r)*sigma(i)/sigma0;
        ePeSum += ePe;
        rSum   += r;
        sigma(i) = 1.;
        if((s > huber) && (r > 1e-4)) // redundancy: it is possible to estimate sigma?
        {
          sigma(i) = std::pow(s/huber, huberPower);
          countOutlierNew++;
        }
      }

      const Double sigma0New = standardDeviation(ePeSum, rSum, huber, huberPower);
      if((countOutlierNew == 0) || ((countOutlier == countOutlierNew) && (std::fabs(sigma0New-sigma0)/sigma0 < 0.001)))
        break;
      sigma0       = sigma0New;
      countOutlier = countOutlierNew;
    } // for(iter)

    return x;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

Matrix Vce::cosTransform(UInt length)
{
  Matrix T(length, length);
  if(!length)
    return T;
  for(UInt i=0; i<length; i++)
    for(UInt k=0; k<length; k++)
      T(i,k) = 2*std::cos(PI*i*k/(length-1));
  T.column(0)        *= 0.5;
  T.column(length-1) *= 0.5;
  T *= 1./std::sqrt(2.*(length-1)); // normalize
  return T;
}

/***********************************************/

Matrix Vce::readCovarianceFunction(const FileName &name, UInt length, UInt columns, Double sampling)
{
  try
  {
    if(!name.empty())
    {
      Matrix covFunc;
      readFileMatrix(name, covFunc);
      if(covFunc.columns() != columns+1)
        throw(Exception("input apriori covariance function <"+name.str()+"> seems not to be compatible"));
      if(covFunc.rows() < length)
        throw(Exception("input apriori covariance function <"+name.str()+"> is too short"));
      if(covFunc.rows() && (std::fabs(covFunc(1,0)-covFunc(0,0)-sampling) > 1e-3)) // test sampling
        throw(Exception("input apriori covariance function <"+name.str()+"> has wrong sampling"));
      return covFunc.row(0,length);
    }

    // default is white noise
    Matrix covFunc(length, 1+columns);
    for(UInt i=0; i<covFunc.rows(); i++)
      covFunc(i,0) = i*sampling;
    if(covFunc.rows())
      for(UInt k=0; k<columns; k++)
        covFunc(0,1+k) = 1.0;
    return covFunc;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Vce::redundancy(const_MatrixSliceRef W, const_MatrixSliceRef We, const_MatrixSliceRef WA, const_MatrixSliceRef WB,
                     Matrix &R, Vector &WWe)
{
  try
  {
    // to compute e^T Sigma^-1 V_i Sigma^-1 e
    // WWe = Sigma^-1 e with Sigma = W^T W
    // (decorrelate again)
    WWe = We;
    triangularSolve(1., W, WWe);

    // Compute redundancy
    // R := Sigma^-1 - Sigma^-1 * (B A) * N^-1 * (B A)^T * Sigma^-1
    //    = Sigma^-1 - W^(-1)*Q1*Q1^T*W^T - W^(-1)*Q2*WA*N^(-1)*N^T*WA^T*Q2^T*W^T
    // ---------------------------------------------------------------------
    R = W;
    cholesky2Inverse(R);
    if(WA.size())
    {
      Matrix WWA = WA;
      triangularSolve(1., W, WWA);
      rankKUpdate(-1., WWA.trans(), R);
    }
    if(WB.size())
    {
      Matrix WWB = WB;
      triangularSolve(1., W, WWB);
      rankKUpdate(-1., WWB.trans(), R);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Vce::matrix(const_MatrixSliceRef R, const_MatrixSliceRef WWe, const_MatrixSliceRef Cov,
                 Double &ePe, Double &redundancy)
{
  try
  {
    Matrix CC = Cov;
    Matrix RR = R;
    fillSymmetric(RR);
    fillSymmetric(CC);
    ePe        += inner(WWe, CC * WWe); // e.T Sigma^-1 C Sigma^-1 e
    redundancy += inner(CC, RR);        // trace(RC)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Vce::psd(const_MatrixSliceRef R, const_MatrixSliceRef WWe,
              const std::vector<UInt> &index, Double sigma, const_MatrixSliceRef CosTransform, const_MatrixSliceRef Psd,
              MatrixSliceRef ePe, MatrixSliceRef redundancy, Double &ePeSum, Double &redundancySum)
{
  try
  {
    const UInt countAxis = Psd.columns();
    for(UInt idAxis=0; idAxis<countAxis; idAxis++)
    {
      Vector e(CosTransform.rows());
      Vector r(CosTransform.rows());
      for(UInt i=0; i<WWe.rows()/countAxis; i++)
        for(UInt k=i; k<WWe.rows()/countAxis; k++)
        {
          e(index.at(k)-index.at(i)) += WWe(countAxis*i+idAxis,0)*WWe(countAxis*k+idAxis,0);
          r(index.at(k)-index.at(i)) += R(countAxis*i+idAxis, countAxis*k+idAxis);
        }
      e.row(1, e.rows()-1) *= 2.; // consider lower triangular of matrix
      r.row(1, r.rows()-1) *= 2.;

      for(UInt idFreq=0; idFreq<Psd.rows(); idFreq++)
      {
        const_MatrixSliceRef cov(CosTransform.column(idFreq));
        const Double ePeTmp        = pow(sigma,2) * Psd(idFreq, idAxis) * inner(e, cov);
        const Double redundancyTmp = pow(sigma,2) * Psd(idFreq, idAxis) * inner(r, cov);
        if((ePeTmp > 0) && !std::isnan(ePeTmp) && (redundancyTmp > 0) && !std::isnan(redundancyTmp))
        {
          ePe(idFreq, idAxis)        += ePeTmp;
          redundancy(idFreq, idAxis) += redundancyTmp;
          ePeSum                     += ePeTmp;
          redundancySum              += redundancyTmp;
        }
      } // for(idFreq)
    } // for(idAxis)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Vce::estimatePsd(MatrixSliceRef ePe, MatrixSliceRef redundancy, MatrixSliceRef Psd, Double &maxFactor, Bool jointZeroFrequency)
{
  try
  {
    if(!Psd.size())
      return;

    // freq=0, special case: mean is usually removed -> not estimable
    // -> estimate first two frequencies together
    if(jointZeroFrequency)
    {
      ePe.row(1)        += ePe.row(0);
      redundancy.row(1) += redundancy.row(0);
      copy(ePe.row(1),        ePe.row(0));
      copy(redundancy.row(1), redundancy.row(0));
    }

    for(UInt idAxis=0; idAxis<Psd.columns(); idAxis++)
      for(UInt idFreq=0; idFreq<Psd.rows(); idFreq++) // frequencies 0, 1,2,3 ...
      {
        Double factor = ePe(idFreq, idAxis)/redundancy(idFreq, idAxis);
        if(std::isnan(factor) || (factor <= 0))
        {
          logWarning<<idFreq<<". frequency, (idAxis="<<idAxis<<") negative factor = "<<ePe(idFreq, idAxis)<<" / "<<redundancy(idFreq, idAxis)<<Log::endl;
          factor = 1;
        }
        Psd(idFreq, idAxis) *= factor;
        maxFactor = std::max(maxFactor, sqrt(exp(fabs(log(factor)))));
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double Vce::meanSigma(const Vector &sigma)
{
  try
  {
    std::vector<Double> s = sigma;
    std::sort(s.begin(), s.end());
    s.erase(std::remove_if(s.begin(), s.end(), [](Double x){return x <= 0;}), s.end());
    if(!s.size())
      return 0;

    const UInt begin  = s.size()/4;
    const UInt end    = s.size()-s.size()/4;
    Double     sigma0 = 0;
    for(UInt i=begin; i<end; i++)
      sigma0 += s.at(i)/(end-begin);
    return sigma0;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
