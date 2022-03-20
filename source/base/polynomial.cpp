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

void Polynomial::init(const std::vector<Time> &times, UInt degree, Bool throwException,
                      Bool leastSquares, Double range, Double extrapolation, Double margin)
{
  try
  {
    this->times          = times;
    this->degree         = degree;
    this->throwException = throwException;
    this->isLeastSquares = leastSquares;
    this->sampling       = medianSampling(times).seconds();
    this->range          = range * ((range < 0) ? -sampling : 1.);
    this->extrapolation  = extrapolation * ((extrapolation < 0) ? -sampling : 1.);

    if(times.size() < degree+1)
      throw(Exception("Not enough data points ("+times.size()%"%i) to interpolate with polynomial degree "s+degree%"%i"s));
    if(std::adjacent_find(times.begin(), times.end(), [margin](const Time &t1, const Time &t2){return (t2-t1).seconds() <= margin;}) != times.end())
      throw(Exception("Input time series is unordered or contains duplicates"));

    std::vector<Bool> isConstInterval(times.size()-1);
    for(UInt i=0; i<isConstInterval.size(); i++)
      isConstInterval.at(i) = (std::fabs((times.at(i+1)-times.at(i)).seconds()-sampling) < margin);

    isPrecomputed.clear();
    isPrecomputed.resize(times.size()-degree, TRUE);
    if(degree > 0)
      for(UInt i=0; i<isPrecomputed.size(); i++)
        isPrecomputed.at(i) = !std::any_of(isConstInterval.begin()+i, isConstInterval.begin()+(i+degree), [](Bool b){return !b;});

    // precomputed polynomial interpolation matrix
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

Matrix Polynomial::interpolate(const std::vector<Time> &timesNew, const_MatrixSliceRef A, UInt rowsPerEpoch, UInt derivative) const
{
  try
  {
    if(!derivative && !isLeastSquares && (timesNew == times)) // need interpolation?
      return A;

    Matrix B(rowsPerEpoch*timesNew.size(), A.columns());
    UInt count = degree+1;
    auto searchStart = times.begin();
    for(UInt i=0; i<timesNew.size(); i++)
    {
      // find interval
      // -------------
      UInt idx;
      if(!isLeastSquares)
      {
        searchStart = std::upper_bound(searchStart, times.end(), timesNew.at(i)); // first epoch greater than interpolation point
        idx = std::min(std::max(static_cast<UInt>(std::distance(times.begin(), searchStart)), count)-count, times.size()-count);

        auto centricity = [&](UInt idx) {return std::max(std::fabs((timesNew.at(i)-times.at(idx)).seconds()),
                                                         std::fabs((timesNew.at(i)-times.at(idx+degree)).seconds()));};
        Double c1 = centricity(idx), c2;
        for(; (idx+count<times.size()) && ((c2 = centricity(idx+1)) < c1); idx++)
          c1 = c2;

        if((times.at(idx+degree)-times.at(idx)).seconds() > range) // polynomial data points not within range
        {
          B.row(i*rowsPerEpoch, rowsPerEpoch).fill(NAN_EXPR);
          if(throwException)
            throw(Exception("cannot interpolate at "+timesNew.at(i).dateTimeStr()));
          continue;
        }
      }
      else
      {
        auto searchEnd = std::upper_bound(searchStart, times.end(), timesNew.at(i)+seconds2time(range)); // first epoch outside search interval
        searchStart    = std::lower_bound(searchStart, searchEnd,   timesNew.at(i)-seconds2time(range)); // first epoch greater or equal than search interval
        count = static_cast<UInt>(std::distance(searchStart, searchEnd));
        if(count < degree+1) // not enough points
        {
          B.row(i*rowsPerEpoch, rowsPerEpoch).fill(NAN_EXPR);
          if(throwException)
            throw(Exception("not enough points to fit at "+timesNew.at(i).dateTimeStr()));
          continue;
        }
        idx = std::distance(times.begin(), searchStart);
      }

      // check if we are allowed to extrapolate
      // --------------------------------------
      if(((times.at(idx) - timesNew.at(i)).seconds()         > extrapolation) || // all points are after newTime and we are too far away
         ((timesNew.at(i) - times.at(idx+count-1)).seconds() > extrapolation))   // all points are before newTime and we are too far away
      {
        B.row(i*rowsPerEpoch, rowsPerEpoch).fill(NAN_EXPR);
        if(throwException)
          throw(Exception("cannot extrapolate at "+timesNew.at(i).dateTimeStr()));
        continue;
      }

      // compute interpolation coefficients
      // ----------------------------------
      Vector coeff;
      if(!isLeastSquares && isPrecomputed.at(idx))
      {
        const Double tau = (timesNew.at(i)-times.at(idx)).seconds()/sampling - degree/2.;
        coeff = Vector(W.rows());
        Double t = 1.0;
        for(UInt n=derivative; n<=degree; n++)
        {
          Double d = 1.0;
          for(UInt k=0; k<derivative; k++)
            d *= (n-k)/sampling;
          axpy(d*t, W.column(n), coeff);
          t *= tau;
        }
      }
      else
      {
        // polynomial matrix
        Matrix P(count, degree+1);
        for(UInt k=0; k<count; k++)
        {
          const Double factor = (timesNew.at(i)-times.at(idx+k)).seconds()/sampling;
          P(k,0) = 1.0;
          for(UInt n=1; n<=degree; n++)
            P(k,n) = factor * P(k,n-1);
        }

        coeff = Vector(count);
        coeff(derivative) = 1.;
        for(UInt n=1; n<=derivative; n++)
          coeff(derivative) *= n/sampling;
        if(P.rows() > P.columns()) // solve with QR-decomposition
        {
          const Vector tau = QR_decomposition(P);
          triangularSolve(1., P.row(0, P.columns()).trans(), coeff.row(0, P.columns())); // R^(-T)*coeff
          QMult(P, tau, coeff);                                                          // coeff := Q*R^(-T)*coeff
        }
        else
        {
          solveInPlace(Matrix(P.trans()), coeff);
        }
      }

      // interpolate
      // -----------
      if(rowsPerEpoch == 1)
        matMult(1., coeff.trans(), A.row(idx, coeff.rows()), B.row(i, rowsPerEpoch));
      else
        for(UInt k=0; k<coeff.rows(); k++)
          axpy(coeff(k), A.row(rowsPerEpoch*(idx+k), rowsPerEpoch), B.row(rowsPerEpoch*i, rowsPerEpoch));
    }

    return B;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
