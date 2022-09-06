/***********************************************/
/**
* @file interpolatorTimeSeriesFillGapsLeastSquaresPolynomialFit.h
*
* @brief Gap filling with least squares polynomial fit.
*
* @author Torsten Mayer-Guerr
* @date 2022-08-01
*
*/
/***********************************************/

#ifndef __GROOPS_INTERPOLATORTIMESERIESFILLGAPSLEASTSQUARESPOLYNOMIALFIT__
#define __GROOPS_INTERPOLATORTIMESERIESFILLGAPSLEASTSQUARESPOLYNOMIALFIT__

// Latex documentation
#ifdef DOCSTRING_InterpolatorTimeSeries
static const char *docstringInterpolatorTimeSeriesFillGapsLeastSquaresPolynomialFit = R"(
\subsection{Fill gaps with least squares polynomial fit}
)";
#endif

/***********************************************/

#include "base/polynomial.h"
#include "classes/interpolatorTimeSeries/interpolatorTimeSeries.h"

/***** CLASS ***********************************/

/** @brief Interpolation using a least squares polynomial fit.
* @ingroup interpolatorTimeSeriesGroup
* @see InterpolatorTimeSeries */
class InterpolatorTimeSeriesFillGapsLeastSquaresPolynomialFit : public InterpolatorTimeSeries
{
  UInt              degree;
  Double            maxDataGap, maxDataSpan;
  Double            margin;
  std::vector<Time> times;
  Bool              throwException;

public:
  InterpolatorTimeSeriesFillGapsLeastSquaresPolynomialFit(Config &config);

  void init(const std::vector<Time> &times, Bool throwException) override;
  Matrix interpolate(const std::vector<Time> &timesNew, const_MatrixSliceRef A, UInt rowsPerEpoch) const override;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline InterpolatorTimeSeriesFillGapsLeastSquaresPolynomialFit::InterpolatorTimeSeriesFillGapsLeastSquaresPolynomialFit(Config &config)
{
  try
  {
    readConfig(config, "polynomialDegree", degree,      Config::DEFAULT, "3",    "degree of the estimated polynomial");
    readConfig(config, "maxDataGap",       maxDataGap,  Config::DEFAULT, "100",  "[seconds] max data gap to interpolate");
    readConfig(config, "maxDataSpan",      maxDataSpan, Config::DEFAULT, "20",   "[seconds] time span on each side used for least squares fit");
    readConfig(config, "margin",           margin,      Config::DEFAULT, "1e-5", "[seconds] margin for identical times");
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void InterpolatorTimeSeriesFillGapsLeastSquaresPolynomialFit::init(const std::vector<Time> &times, Bool throwException)
{
  try
  {
    this->times = times;
    this->throwException = throwException;

    const Double sampling = medianSampling(times).seconds();
    if(maxDataGap  < 0) maxDataGap  *= -sampling;
    if(maxDataSpan < 0) maxDataSpan *= -sampling;
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Matrix InterpolatorTimeSeriesFillGapsLeastSquaresPolynomialFit::interpolate(const std::vector<Time> &timesNew, const_MatrixSliceRef A, UInt rowsPerEpoch) const
{
  try
  {
    if(timesNew == times) // need interpolation?
      return A;
    Matrix B(rowsPerEpoch*timesNew.size(), A.columns());

    // new points before data?
    // -----------------------
    UInt iStart=0;
    for(; (iStart<timesNew.size()) && ((timesNew.at(iStart)-times.front()).seconds() < -margin); iStart++)
    {
      if(throwException)
        throw(Exception("extrapolation at "+timesNew.front().dateTimeStr()+" before start of data at "+times.front().dateTimeStr()));
      B.row(iStart*rowsPerEpoch, rowsPerEpoch).fill(NAN_EXPR);
    }

    // new points after data?
    // ----------------------
    UInt iEnd = timesNew.size();
    for(; (iEnd-- > 0) && ((timesNew.at(iEnd)-times.back()).seconds() > margin);)
    {
      if(throwException)
        throw(Exception("extrapolation at "+timesNew.back().dateTimeStr()+" after end of data at "+times.back().dateTimeStr()));
      B.row(iEnd*rowsPerEpoch, rowsPerEpoch).fill(NAN_EXPR);
    }

    UInt k = 0;
    for(UInt i=iStart; i<=iEnd; i++)
    {
      for(; (k+1 < times.size()) && ((times.at(k+1)-timesNew.at(i)).seconds() <= margin); k++);

      // no gap?
      if(std::fabs((timesNew.at(i)-times.at(k)).seconds()) <= margin)
      {
        copy(A.row(k, rowsPerEpoch), B.row(i, rowsPerEpoch));
        continue;
      }

      if((times.at(k+1)-times.at(k)).seconds() > maxDataGap)
      {
        if(throwException)
          throw(Exception("gap too large between "+times.at(k).dateTimeStr()+" - "+times.at(k+1).dateTimeStr()));
        B.row(i*rowsPerEpoch, rowsPerEpoch).fill(NAN_EXPR);
        continue;
      }

      // find interval
      // -------------
      UInt idx = k;         // left interval
      for(; (idx > 0) && ((times.at(k)-times.at(idx-1)).seconds() <= maxDataSpan); idx--);
      UInt count = k-idx+1; // right interval
      for(; (idx+count < times.size()) && ((times.at(idx+count)-times.at(k+1)).seconds() <= maxDataSpan); count++);
      if(count < degree+1)  // not enough points
      {
        if(throwException)
          throw(Exception("not enough points to fit at "+timesNew.at(i).dateTimeStr()));
        B.row(i*rowsPerEpoch, rowsPerEpoch).fill(NAN_EXPR);
        continue;
      }

      // compute interpolation coefficients
      // ----------------------------------
      Matrix P(count, degree+1); // polynomial matrix
      for(UInt k=0; k<count; k++)
      {
        const Double factor = (timesNew.at(i)-times.at(idx+k)).seconds()/(times.at(idx+count-1)-times.at(idx)).seconds();
        P(k,0) = 1.0;
        for(UInt n=1; n<=degree; n++)
          P(k,n) = factor * P(k,n-1);
      }

      Vector coeff(count);
      coeff(0) = 1.;
      const Vector tau = QR_decomposition(P);
      triangularSolve(1., P.row(0, P.columns()).trans(), coeff.row(0, P.columns())); // R^(-T)*coeff
      QMult(P, tau, coeff);                                                          // coeff := Q*R^(-T)*coeff

      // interpolate
      // -----------
      for(UInt k=0; k<coeff.rows(); k++)
        axpy(coeff(k), A.row(rowsPerEpoch*(idx+k), rowsPerEpoch), B.row(rowsPerEpoch*i, rowsPerEpoch));
    } // for(i)

    return B;
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
