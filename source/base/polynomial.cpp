/***********************************************/
/**
* @file polynomial.cpp
*
* @brief Interpolation by polynomial.
*
* @author Torsten Mayer-Guerr
* @date 2017-05-27
*
*/
/***********************************************/

#include "base/importStd.h"
#include "base/polynomial.h"

/***********************************************/

void Polynomial::init(UInt degree)
{
  try
  {
    if(W.rows() == degree+1)
      return;

    // polynomial interpolation matrix
    W = Matrix(degree+1, degree+1);
    for(UInt i=0; i<W.rows(); i++)
    {
      W(0,i) = 1.0;
      for(UInt n=1; n<W.columns(); n++)
        W(n,i) = (i-degree/2.) * W(n-1,i);
    }
    inverse(W);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix Polynomial::interpolate(const std::vector<Time> &timesNew, const std::vector<Time> &times, const_MatrixSliceRef A, UInt rowsPerEpoch) const
{
  try
  {
    if(timesNew == times) // need interpolation?
      return A;

    if(!W.size())
      throw(Exception("Polynomial is not initialized"));
    if(times.size()*rowsPerEpoch != A.rows())
      throw(Exception("Dimension error: time.size()="+times.size()%"%i and A("s+A.rows()%"%i x "s+A.columns()%"%i)"s));
    if(times.size()<W.rows())
      throw(Exception("Too few epochs with A("+A.rows()%"%i x "s+A.columns()%"%i)"s));

    const UInt   degree = W.rows()-1;
    const Time   time0  = times.at(0);
    const Double dt     = (times.at(1) - times.at(0)).seconds(); // assume constant sampling

    Matrix B(rowsPerEpoch*timesNew.size(), A.columns());
    for(UInt i=0; i<timesNew.size(); i++)
    {
      // interpolation interval
      const UInt   idx = std::min(static_cast<UInt>(std::max(round((timesNew.at(i)-time0).seconds()/dt-degree/2.), 0.)), times.size()-degree-1);
      const Double tau = (timesNew.at(i)-times.at(idx)).seconds()/dt - degree/2.;

      // interpolation coefficients
      Vector coeff(degree+1);
      Double factor = 1.0;
      for(UInt n=0; n<coeff.rows(); n++)
      {
        axpy(factor, W.column(n), coeff);
        factor *= tau;
      }

      // interpolate
      MatrixSlice Sum(B.row(rowsPerEpoch*i, rowsPerEpoch));
      if(rowsPerEpoch==1)
        matMult(1., coeff.trans(), A.row(rowsPerEpoch*idx, coeff.rows()), Sum);
      else
        for(UInt k=0; k<coeff.rows(); k++)
          axpy(coeff(k), A.row(rowsPerEpoch*(idx+k), rowsPerEpoch), Sum);
    }

    return B;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix Polynomial::derivative(Double sampling, const_MatrixSliceRef A, UInt rowsPerEpoch) const
{
  try
  {
    if(!W.size())
      throw(Exception("Polynomial is not initialized"));
    if(A.rows()<rowsPerEpoch*W.rows())
      throw(Exception("Too few epochs with A("+A.rows()%"%i x "s+A.columns()%"%i)"s));

    const UInt degree = W.rows()-1;
    const UInt size   = A.rows()/rowsPerEpoch;
    Vector coeff;
    Matrix B(A.rows(), A.columns());
    for(UInt i=0; i<size; i++)
    {
      const UInt idx = std::min(std::max(i, degree/2)-degree/2, size-degree-1); // interval
      if((i<=degree/2) || (i>=size-(degree+1)/2))
      {
        coeff = Vector(degree+1);  // derivation coefficients
        const Double tau = i-idx-degree/2.;
        Double factor = 1./sampling;
        for(UInt n=1; n<coeff.rows(); n++)
        {
          axpy(n*factor, W.column(n), coeff);
          factor *= tau;
        }
      }

      // interpolate
      MatrixSlice Sum(B.row(rowsPerEpoch*i, rowsPerEpoch));
      if(rowsPerEpoch==1)
        matMult(1., coeff.trans(), A.row(rowsPerEpoch*idx, coeff.rows()), Sum);
      else
        for(UInt k=0; k<coeff.rows(); k++)
          axpy(coeff(k), A.row(rowsPerEpoch*(idx+k), rowsPerEpoch), Sum);
    }

    return B;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix Polynomial::derivative2nd(Double sampling, const_MatrixSliceRef A, UInt rowsPerEpoch) const
{
  try
  {
    if(!W.size())
      throw(Exception("Polynomial is not initialized"));
    if(A.rows()<rowsPerEpoch*W.rows())
      throw(Exception("Too few epochs with A("+A.rows()%"%i x "s+A.columns()%"%i)"s));

    const UInt degree = W.rows()-1;
    const UInt size   = A.rows()/rowsPerEpoch;
    Vector coeff;
    Matrix B(A.rows(), A.columns());
    for(UInt i=0; i<size; i++)
    {
      const UInt idx = std::min(std::max(i, degree/2)-degree/2, size-degree-1); // interval
      if((i<=degree/2) || (i>=size-(degree+1)/2))
      {
        coeff = Vector(degree+1);  // derivation coefficients
        const Double tau = i-idx-degree/2.;
        Double factor = 1./sampling/sampling;
        for(UInt n=2; n<coeff.rows(); n++)
        {
          axpy(n*(n-1)*factor, W.column(n), coeff);
          factor *= tau;
        }
      }

      // interpolate
      MatrixSlice Sum(B.row(rowsPerEpoch*i, rowsPerEpoch));
      if(rowsPerEpoch==1)
        matMult(1., coeff.trans(), A.row(rowsPerEpoch*idx, coeff.rows()), Sum);
      else
        for(UInt k=0; k<coeff.rows(); k++)
          axpy(coeff(k), A.row(rowsPerEpoch*(idx+k), rowsPerEpoch), Sum);
    }

    return B;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix Polynomial::integration(Double sampling, const_MatrixSliceRef A, UInt rowsPerEpoch) const
{
  try
  {
    if(!W.size())
      throw(Exception("Polynomial is not initialized"));
    if(A.rows()<rowsPerEpoch*W.rows())
      throw(Exception("Too few epochs with A("+A.rows()%"%i x "s+A.columns()%"%i)"s));

    const UInt degree = W.rows()-1;
    const UInt size   = A.rows()/rowsPerEpoch;
    Vector coeff;
    Matrix B(A.rows(), A.columns());
    for(UInt i=1; i<size; i++)
    {
      const UInt idx = std::min(std::max(i, degree/2)-degree/2, size-degree-1); // interval
      if((i<=degree/2) || (i>=size-(degree+1)/2))
      {
        coeff = Vector(degree+1);  // derivation coefficients
        const Double tau1 = i-idx-degree/2.-1;
        const Double tau2 = i-idx-degree/2.;
        for(UInt n=0; n<coeff.rows(); n++)
          axpy(sampling/(n+1)*(std::pow(tau2, n+1)-std::pow(tau1, n+1)), W.column(n), coeff);
      }

      // interpolate
      MatrixSlice Sum(B.row(rowsPerEpoch*i, rowsPerEpoch));
      Sum += B.row(rowsPerEpoch*(i-1), rowsPerEpoch);
      if(rowsPerEpoch==1)
        matMult(1., coeff.trans(), A.row(rowsPerEpoch*idx, coeff.rows()), Sum);
      else
        for(UInt k=0; k<coeff.rows(); k++)
          axpy(coeff(k), A.row(rowsPerEpoch*(idx+k), rowsPerEpoch), Sum);
    }

    return B;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
